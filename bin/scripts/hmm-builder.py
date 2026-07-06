#!/usr/bin/env python
"""
Build per-cluster HMMs from a protein fasta + cluster membership table.

Cutoff calibration methods:
  - 'simple' (default): NC = max bit score on non-members, TC = min bit
               score on members. GA placed between NC and TC with an
               epsilon margin. Fast — one hmmsearch per cluster.
  - 'kfold'  : k-fold CV with F-measure maximization. For each cluster,
               split positives into k folds; for each fold, build an HMM
               from the other folds, score held-out positives (TP) vs all
               non-member proteins (FP), find T maximizing F1. GA = mean
               of per-fold T's. Much more expensive (k rebuilds + k
               hmmsearches per cluster).

Notes:
  - Singleton clusters (n=1) are excluded by default. Use --include_singletons
    to build them as single-sequence (phmmer-equivalent) HMMs. Single-sequence
    HMMs cannot capture family-level conservation and their calibration is
    degenerate (TC = self-score), so they are not best practice for an HMM
    database. When included, they are always calibrated with the 'simple'
    method regardless of --calibration_method.
  - Default output is stdout; logs go to stderr (loguru).
"""
import argparse
import gc
import gzip
import hashlib
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyfastx
from loguru import logger
from tqdm import tqdm
from pyfamsa import Aligner, Sequence
from pytrimal import Alignment, AutomaticTrimmer
import pyhmmer
from pyhmmer.easel import (
    Alphabet,
    TextSequence,
    TextMSA,
)
from pyhmmer.plan7 import (
    Builder,
    Background,
    Pipeline,
)


# ---------------------------------------------------------------------------
# IO
# ---------------------------------------------------------------------------

# Canonical 20 + X (unknown); anything else gets replaced with X at load time
# to avoid crashes in pyfamsa (MIQS is 20x20) and pytrimal (similarity matrices
# throw std::out_of_range on B/Z/J/O/U). pyhmmer handles degenerates correctly
# but we sanitize upstream so the whole pipeline sees a consistent alphabet.
_CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWYX")
_AA_TRANSLATE_TABLE = str.maketrans(
    {ch: "X" for ch in "BZJOU*bzjou"}
)


def sanitize_sequence(seq):
    """Uppercase and replace non-canonical residues with X.

    Anything not in the 20 canonical + X becomes X. This includes
    B (Asx), Z (Glx), J (Xle), O (pyrrolysine), U (selenocysteine),
    '*' (stop codon), and stray lowercase / whitespace.
    """
    seq_upper = seq.upper().translate(_AA_TRANSLATE_TABLE)
    # Catch anything else that isn't canonical (whitespace, digits, punctuation)
    return "".join(ch if ch in _CANONICAL_AA else "X" for ch in seq_upper)


def load_proteins(fasta_path):
    """Return dict id_protein -> sanitized sequence string."""
    logger.info(f"Loading proteins from {fasta_path}")
    protein_to_seq = {}
    n_sanitized = 0
    iterator = tqdm(
        pyfastx.Fasta(fasta_path, build_index=False),
        desc="Loading proteins",
        unit="seq",
        unit_scale=True,
        mininterval=1.0,
    )
    for name, seq in iterator:
        pid = name.split()[0]
        sanitized = sanitize_sequence(seq)
        if sanitized != seq:
            n_sanitized += 1
        protein_to_seq[pid] = sanitized
    logger.info(f"Loaded {len(protein_to_seq)} proteins ({n_sanitized} had non-canonical residues replaced with X)")
    return protein_to_seq


def load_clusters(clusters_path):
    """Return dict id_protein_cluster -> list of id_protein.

    Uses groupby rather than iterrows, which is orders of magnitude faster
    for large tables (6.9M rows drops from ~8 min to a few seconds).
    """
    logger.info(f"Loading cluster table from {clusters_path}")
    df_clusters = pd.read_csv(
        clusters_path,
        sep="\t",
        header=None,
        names=["id_protein", "id_protein_cluster"],
        dtype=str,
    )
    logger.info(f"Grouping {len(df_clusters)} assignments by cluster")
    groups = df_clusters.groupby("id_protein_cluster", sort=False)["id_protein"]
    cluster_to_proteins = {
        cluster_id: proteins.tolist()
        for cluster_id, proteins in tqdm(
            groups,
            desc="Grouping clusters",
            total=groups.ngroups,
            unit="cluster",
            mininterval=1.0,
        )
    }
    logger.info(
        f"Loaded {len(cluster_to_proteins)} clusters covering {len(df_clusters)} protein-cluster assignments"
    )
    return cluster_to_proteins


# ---------------------------------------------------------------------------
# MSA / trim / build
# ---------------------------------------------------------------------------

def align_sequences(protein_ids, protein_to_seq, threads):
    """pyfamsa MSA. Returns list of (id, aligned_seq) or None if <2 usable seqs.

    pyfamsa's Aligner has no seed parameter. Default single-linkage guide tree
    is deterministic given input order (Prim's MST, ties broken by input
    order), so reproducibility comes from stable input order, not seeding.
    Medoid-tree heuristic is stochastic and not exposed for seeding — avoid it.
    """
    sequences = [
        Sequence(pid.encode(), protein_to_seq[pid].encode())
        for pid in protein_ids
        if pid in protein_to_seq
    ]
    if len(sequences) < 2:
        return None
    aligner = Aligner(threads=threads)
    msa = aligner.align(sequences)
    return [(seq.id.decode(), seq.sequence.decode()) for seq in msa]


def trim_alignment(aligned_pairs, method):
    """pytrimal. Returns list of (id, trimmed_seq)."""
    names = [pid.encode() for pid, _ in aligned_pairs]
    seqs = [s for _, s in aligned_pairs]
    alignment = Alignment(names, seqs)
    trimmer = AutomaticTrimmer(method=method)
    trimmed = trimmer.trim(alignment)
    return [
        (name.decode(), seq)
        for name, seq in zip(trimmed.names, trimmed.sequences)
    ]


def build_hmm_from_msa(cluster_id, aligned_pairs, alphabet, builder, background):
    """Build HMM from MSA."""
    text_msa = TextMSA(
        name=cluster_id.encode(),
        sequences=[
            TextSequence(name=pid.encode(), sequence=seq)
            for pid, seq in aligned_pairs
        ],
    )
    digital_msa = text_msa.digitize(alphabet)
    hmm, _, _ = builder.build_msa(digital_msa, background)
    return hmm


def build_hmm_from_single_sequence(cluster_id, protein_id, seq_str, alphabet, builder, background):
    """Build a single-sequence HMM (phmmer-equivalent).

    Used for singleton clusters where MSA is undefined. The resulting HMM
    has position-specific emissions derived from the single observed residue
    plus the Dirichlet prior, so it behaves like a PSSM with HMMER's
    statistical framework on top. It will detect close homologs of the
    source sequence but lacks family-level conservation signal.
    """
    text_seq = TextSequence(name=protein_id.encode(), sequence=seq_str)
    digital_seq = text_seq.digitize(alphabet)
    hmm, _, _ = builder.build(digital_seq, background)
    hmm.name = cluster_id.encode()
    return hmm


# Filesystem-safe cluster id for use in filenames
_UNSAFE_FILENAME_CHARS = re.compile(r"[^A-Za-z0-9._-]+")


def sanitize_for_filename(cluster_id):
    """Replace characters unsafe in filenames with underscores."""
    return _UNSAFE_FILENAME_CHARS.sub("_", cluster_id)


def write_msa_gzip(cluster_id, aligned_pairs, msa_directory):
    """Write an aligned FASTA to {msa_directory}/{sanitized_cluster_id}.msa.fasta.gz."""
    safe = sanitize_for_filename(cluster_id)
    path_output = Path(msa_directory) / f"{safe}.msa.fasta.gz"
    with gzip.open(path_output, "wt") as file_handle:
        for pid, seq in aligned_pairs:
            print(f">{pid}", file=file_handle)
            print(seq, file=file_handle)


def make_hmm_for_cluster(
    cluster_id,
    protein_ids,
    protein_to_seq,
    do_trim,
    trim_method,
    min_seqs_for_trim,
    alphabet,
    builder,
    background,
    threads,
    msa_directory=None,
):
    """Full MSA -> (trim) -> HMM pipeline for one cluster. Returns HMM or None.

    If msa_directory is given, writes the final (post-trim if trimmed) MSA
    as gzipped FASTA.
    """
    aligned = align_sequences(protein_ids, protein_to_seq, threads)
    if aligned is None:
        return None
    if do_trim and len(aligned) >= min_seqs_for_trim:
        try:
            aligned = trim_alignment(aligned, trim_method)
        except Exception as e:
            logger.warning(f"Cluster {cluster_id}: trimming failed ({e}), using untrimmed MSA")
    if msa_directory is not None:
        try:
            write_msa_gzip(cluster_id, aligned, msa_directory)
        except Exception as e:
            logger.warning(f"Cluster {cluster_id}: MSA write failed ({e})")
    try:
        return build_hmm_from_msa(cluster_id, aligned, alphabet, builder, background)
    except Exception as e:
        logger.warning(f"Cluster {cluster_id}: HMM build failed ({e})")
        return None


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def digitize_proteome(protein_to_seq, alphabet):
    """Return list of DigitalSequence for all proteins."""
    digital_seqs = []
    for pid, seq in tqdm(
        protein_to_seq.items(),
        desc="Digitizing proteome",
        unit="seq",
        unit_scale=True,
        mininterval=1.0,
    ):
        ts = TextSequence(name=pid.encode(), sequence=seq)
        digital_seqs.append(ts.digitize(alphabet))
    return digital_seqs


def score_hmm_vs_proteome(
    hmm,
    digital_proteome,
    threads,
    force_members=None,
    alphabet=None,
    background=None,
):
    """
    Run hmmsearch of one HMM against digital proteome. Returns dict
    id_protein -> bit score (full-sequence).

    If force_members is given, members missing from default hits are
    force-scored via a permissive Pipeline and merged in.
    """
    hits_iter = pyhmmer.hmmsearch([hmm], digital_proteome, cpus=threads)
    top_hits = next(iter(hits_iter))
    score_dict = {hit.name.decode(): hit.score for hit in top_hits}

    if force_members:
        missing = [m for m in force_members if m not in score_dict]
        if missing:
            missing_set = set(missing)
            subset = [ds for ds in digital_proteome if ds.name.decode() in missing_set]
            if subset:
                try:
                    pipeline = Pipeline(
                        alphabet=alphabet,
                        background=background,
                        E=1e9,
                        domE=1e9,
                        incE=1e9,
                        incdomE=1e9,
                        bias_filter=False,
                        F1=1.0, F2=1.0, F3=1.0,
                    )
                    forced = pipeline.search_hmm(hmm, subset)
                    for hit in forced:
                        name = hit.name.decode()
                        if name not in score_dict:
                            score_dict[name] = hit.score
                except Exception as e:
                    logger.debug(f"Force-scoring failed for HMM {hmm.name.decode()}: {e}")
    return score_dict


# ---------------------------------------------------------------------------
# Calibration: simple
# ---------------------------------------------------------------------------

def calibrate_simple(
    cluster_id,
    hmm,
    members_set,
    digital_proteome,
    threads,
    epsilon,
    force_score_members,
    alphabet,
    background,
):
    """NC = max FP bit score, TC = min TP bit score, GA between them."""
    score_dict = score_hmm_vs_proteome(
        hmm,
        digital_proteome,
        threads,
        force_members=members_set if force_score_members else None,
        alphabet=alphabet,
        background=background,
    )
    tp_scores = [s for name, s in score_dict.items() if name in members_set]
    fp_scores = [s for name, s in score_dict.items() if name not in members_set]

    nc = max(fp_scores) if fp_scores else None
    tc = min(tp_scores) if tp_scores else None

    if nc is None and tc is None:
        return dict(ga=0.0, tc=0.0, nc=0.0, n_tp=0, n_fp=0, note="no_hits")
    if nc is None:
        ga = max(tc - epsilon, 0.0)
        return dict(ga=ga, tc=tc, nc=0.0, n_tp=len(tp_scores), n_fp=0, note="clean")
    if tc is None:
        return dict(ga=nc + epsilon, tc=nc + epsilon, nc=nc,
                    n_tp=0, n_fp=len(fp_scores), note="no_self_hits")
    if tc <= nc:
        return dict(ga=nc + epsilon, tc=tc, nc=nc,
                    n_tp=len(tp_scores), n_fp=len(fp_scores), note="overlap")
    ga = min(nc + epsilon, tc)
    return dict(ga=ga, tc=tc, nc=nc,
                n_tp=len(tp_scores), n_fp=len(fp_scores), note="ok")


# ---------------------------------------------------------------------------
# Calibration: KOfam-style CV
# ---------------------------------------------------------------------------

def fmeasure_optimal_threshold(pos_scores, neg_scores):
    """
    Find bit score T maximizing F1 where predictions are 'positive' for score >= T.
    Returns (T_opt, F_opt, n_pos, n_neg).
    """
    n_pos = len(pos_scores)
    n_neg = len(neg_scores)
    if n_pos == 0:
        return 0.0, 0.0, n_pos, n_neg

    pos = np.asarray(pos_scores, dtype=float)
    neg = np.asarray(neg_scores, dtype=float) if n_neg else np.array([], dtype=float)

    all_scores = np.concatenate([pos, neg]) if n_neg else pos
    candidates = np.unique(np.concatenate([all_scores, [pos.max() + 1.0]]))

    best_t = float(candidates[0])
    best_f = -1.0
    for t in candidates:
        tp = int((pos >= t).sum())
        fp = int((neg >= t).sum()) if n_neg else 0
        fn = n_pos - tp
        if tp == 0:
            f = 0.0
        else:
            precision = tp / (tp + fp)
            recall = tp / (tp + fn)
            f = 2 * precision * recall / (precision + recall)
        if f > best_f:
            best_f = f
            best_t = float(t)
    return best_t, float(best_f), n_pos, n_neg


def _derive_cluster_seed(cluster_id, random_state):
    """Deterministic per-cluster seed derived from global random_state + cluster_id.

    Uses md5 so it is reproducible across Python sessions (unlike hash()).
    usedforsecurity=False lets this work identically on FIPS-mode systems
    where the default md5 constructor is disallowed. Added in Python 3.9.
    """
    h = hashlib.md5(
        f"{random_state}::{cluster_id}".encode(),
        usedforsecurity=False,
    ).hexdigest()
    return int(h[:8], 16)


def calibrate_kfold_cv(
    cluster_id,
    protein_ids,
    protein_to_seq,
    digital_proteome,
    alphabet,
    builder,
    background,
    do_trim,
    trim_method,
    min_seqs_for_trim,
    k_folds,
    threads,
    force_score_members,
    random_state,
):
    """
    k-fold CV calibration. Returns dict with ga/tc/nc + cv info.
    """
    members_set = set(protein_ids)
    n = len(protein_ids)

    effective_k = min(k_folds, n)
    if effective_k < 2:
        return dict(ga=0.0, tc=0.0, nc=0.0, n_tp=0, n_fp=0,
                    note="too_few_for_cv", cv_thresholds=[], cv_fmeasure=[])

    cluster_seed = _derive_cluster_seed(cluster_id, random_state)
    rng = np.random.default_rng(seed=cluster_seed)
    shuffled = list(protein_ids)
    rng.shuffle(shuffled)
    folds = [shuffled[i::effective_k] for i in range(effective_k)]

    cv_thresholds = []
    cv_fmeasures = []
    tc_candidates = []
    nc_candidates = []
    total_tp = 0
    total_fp = 0

    for i in range(effective_k):
        held_out = folds[i]
        training = [p for j, fold in enumerate(folds) if j != i for p in fold]
        if len(training) < 2:
            continue
        fold_hmm = make_hmm_for_cluster(
            cluster_id=f"{cluster_id}__fold{i}",
            protein_ids=training,
            protein_to_seq=protein_to_seq,
            do_trim=do_trim,
            trim_method=trim_method,
            min_seqs_for_trim=min_seqs_for_trim,
            alphabet=alphabet,
            builder=builder,
            background=background,
            threads=threads,
            msa_directory=None,
        )
        if fold_hmm is None:
            continue

        score_dict = score_hmm_vs_proteome(
            fold_hmm,
            digital_proteome,
            threads,
            force_members=set(held_out) if force_score_members else None,
            alphabet=alphabet,
            background=background,
        )
        pos_scores = [score_dict.get(p, 0.0) for p in held_out]
        neg_scores = [s for name, s in score_dict.items() if name not in members_set]

        t_opt, f_opt, n_pos_fold, n_neg_fold = fmeasure_optimal_threshold(pos_scores, neg_scores)
        cv_thresholds.append(t_opt)
        cv_fmeasures.append(f_opt)
        total_tp += n_pos_fold
        total_fp += n_neg_fold
        if pos_scores:
            tc_candidates.append(min(pos_scores))
        if neg_scores:
            nc_candidates.append(max(neg_scores))

    if not cv_thresholds:
        return dict(ga=0.0, tc=0.0, nc=0.0, n_tp=0, n_fp=0,
                    note="cv_failed", cv_thresholds=[], cv_fmeasure=[])

    ga = float(np.mean(cv_thresholds))
    tc = float(np.mean(tc_candidates)) if tc_candidates else ga
    nc = float(np.mean(nc_candidates)) if nc_candidates else 0.0

    note = "ok"
    if tc < ga:
        note = "tc_lt_ga"
    if nc > ga:
        note = "nc_gt_ga"

    return dict(
        ga=ga, tc=tc, nc=nc,
        n_tp=total_tp, n_fp=total_fp,
        note=note,
        cv_thresholds=cv_thresholds,
        cv_fmeasure=cv_fmeasures,
    )


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def build_all_hmms(
    cluster_to_proteins,
    protein_to_seq,
    do_trim,
    trim_method,
    min_seqs_for_trim,
    alphabet,
    builder,
    background,
    threads,
    include_singletons=False,
    msa_directory=None,
    progress_file=None,
):
    """Build one HMM per cluster.

    Singletons (n=1) are excluded by default because single-sequence HMMs
    cannot capture family-level conservation and behave essentially like
    phmmer against the source sequence — they don't add information over
    BLAST/phmmer while inflating the multiple-testing burden. Set
    include_singletons=True to build them anyway as single-sequence HMMs
    via pyhmmer's Builder.build().

    If progress_file is given, the current cluster_id is written to it
    (overwriting) at the start of each iteration. This is a crash-
    diagnostic aid: if a C++ binding (pyfamsa / pytrimal / pyhmmer)
    calls std::terminate() on malformed input, Python's try/except
    cannot catch it — the process is SIGABRT'd mid-iteration. The
    progress file preserves the last cluster attempted so the offender
    is identifiable without rerunning the pipeline.

    Returns (cluster_to_hmm, singleton_cluster_ids). singleton_cluster_ids
    is empty when include_singletons=False.
    """
    cluster_to_hmm = {}
    singleton_cluster_ids = set()
    n_singletons_included = 0
    n_singletons_excluded = 0
    for cluster_id, protein_ids in tqdm(
        cluster_to_proteins.items(),
        desc="Building HMMs",
        total=len(cluster_to_proteins),
        unit="cluster",
        mininterval=1.0,
    ):
        if progress_file is not None:
            # Overwrite each iteration; on crash this is the last cluster
            # attempted. open+write+close to guarantee an fsync-free flush.
            try:
                Path(progress_file).write_text(f"{cluster_id}\n")
            except OSError:
                pass
        present = [p for p in protein_ids if p in protein_to_seq]
        if len(present) == 0:
            logger.warning(f"Cluster {cluster_id}: no members present in fasta — skipping")
            continue
        if len(present) == 1:
            if not include_singletons:
                n_singletons_excluded += 1
                continue
            # Build single-sequence HMM
            pid = present[0]
            try:
                hmm = build_hmm_from_single_sequence(
                    cluster_id, pid, protein_to_seq[pid], alphabet, builder, background,
                )
                cluster_to_hmm[cluster_id] = hmm
                singleton_cluster_ids.add(cluster_id)
                n_singletons_included += 1
                # Write "MSA" (single seq) if requested, for consistency
                if msa_directory is not None:
                    try:
                        write_msa_gzip(cluster_id, [(pid, protein_to_seq[pid])], msa_directory)
                    except Exception as e:
                        logger.warning(f"Cluster {cluster_id}: MSA write failed ({e})")
            except Exception as e:
                logger.warning(f"Cluster {cluster_id}: singleton HMM build failed ({e})")
            continue
        hmm = make_hmm_for_cluster(
            cluster_id=cluster_id,
            protein_ids=present,
            protein_to_seq=protein_to_seq,
            do_trim=do_trim,
            trim_method=trim_method,
            min_seqs_for_trim=min_seqs_for_trim,
            alphabet=alphabet,
            builder=builder,
            background=background,
            threads=threads,
            msa_directory=msa_directory,
        )
        if hmm is not None:
            cluster_to_hmm[cluster_id] = hmm
    if include_singletons:
        logger.info(
            f"Built {len(cluster_to_hmm)} HMMs "
            f"({n_singletons_included} from singletons as single-sequence HMMs)"
        )
    else:
        logger.info(f"Built {len(cluster_to_hmm)} HMMs; excluded {n_singletons_excluded} singleton(s)")
    return cluster_to_hmm, singleton_cluster_ids


def calibrate_all(
    method,
    cluster_to_hmm,
    cluster_to_proteins,
    protein_to_seq,
    digital_proteome,
    alphabet,
    builder,
    background,
    do_trim,
    trim_method,
    min_seqs_for_trim,
    k_folds,
    epsilon,
    force_score_members,
    threads,
    random_state,
    singleton_cluster_ids=None,
):
    """Per-cluster calibration. Attaches cutoffs. Returns df_stats.

    Singletons always use the 'simple' method because kfold CV is undefined
    for n=1. For a singleton: TC = self-score, NC = max non-member score,
    GA between them.
    """
    if singleton_cluster_ids is None:
        singleton_cluster_ids = set()
    cluster_to_members = {c: set(p) for c, p in cluster_to_proteins.items()}
    rows = []

    for cluster_id, hmm in tqdm(
        cluster_to_hmm.items(),
        desc=f"Calibrating cutoffs ({method})",
        total=len(cluster_to_hmm),
        unit="hmm",
        mininterval=1.0,
    ):
        members = cluster_to_members[cluster_id]
        is_singleton = cluster_id in singleton_cluster_ids
        effective_method = "simple" if is_singleton else method
        if effective_method == "simple":
            result = calibrate_simple(
                cluster_id=cluster_id,
                hmm=hmm,
                members_set=members,
                digital_proteome=digital_proteome,
                threads=threads,
                epsilon=epsilon,
                force_score_members=force_score_members,
                alphabet=alphabet,
                background=background,
            )
            if is_singleton:
                # Overwrite the note to make singleton status explicit in stats
                result["note"] = f"singleton_{result['note']}"
        elif effective_method == "kfold":
            result = calibrate_kfold_cv(
                cluster_id=cluster_id,
                protein_ids=[p for p in cluster_to_proteins[cluster_id] if p in protein_to_seq],
                protein_to_seq=protein_to_seq,
                digital_proteome=digital_proteome,
                alphabet=alphabet,
                builder=builder,
                background=background,
                do_trim=do_trim,
                trim_method=trim_method,
                min_seqs_for_trim=min_seqs_for_trim,
                k_folds=k_folds,
                threads=threads,
                force_score_members=force_score_members,
                random_state=random_state,
            )
        else:
            raise ValueError(f"Unknown calibration method: {method}")

        ga, tc, nc = result["ga"], result["tc"], result["nc"]
        hmm.cutoffs.gathering = (ga, ga)
        hmm.cutoffs.trusted = (tc, tc)
        hmm.cutoffs.noise = (nc, nc)

        row = {
            "id_protein_cluster": cluster_id,
            "n_members": len(members),
            "ga": ga,
            "tc": tc,
            "nc": nc,
            "n_tp_hits": result["n_tp"],
            "n_fp_hits": result["n_fp"],
            "note": result["note"],
        }
        if effective_method == "kfold":
            row["cv_thresholds"] = ",".join(f"{t:.3f}" for t in result.get("cv_thresholds", []))
            row["cv_fmeasure"] = ",".join(f"{f:.3f}" for f in result.get("cv_fmeasure", []))
        rows.append(row)

    df_stats = pd.DataFrame(rows)
    return df_stats


def write_hmms(cluster_to_hmm, output):
    """Concatenate all HMMs. `output` is either '-' for stdout or a file path."""
    if output == "-":
        f = sys.stdout.buffer
        close = False
    else:
        f = open(output, "wb")
        close = True
    try:
        for hmm in cluster_to_hmm.values():
            hmm.write(f)
    finally:
        if close:
            f.close()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build per-cluster HMMs from a protein fasta and cluster table, "
            "with optional KOfam-style CV cutoff calibration."
        ),
    )
    parser.add_argument("-f", "--fasta", required=True, help="Protein fasta file")
    parser.add_argument("-c", "--clusters", required=True,
                        help="TSV (no header): id_protein\\tid_protein_cluster")
    parser.add_argument("-o", "--output_hmm", default="-",
                        help="Output concatenated HMM file (default: stdout)")
    parser.add_argument("--output_stats", default=None,
                        help="Optional TSV with per-cluster cutoff stats")
    parser.add_argument("--msa_directory", default=None,
                        help="Optional directory to write per-cluster MSAs "
                             "as gzipped FASTA ({id_protein_cluster}.msa.fasta.gz). "
                             "Written post-trim if trimming is enabled.")
    parser.add_argument("--progress_file", default=None,
                        help="Optional file to record the cluster currently being "
                             "processed (overwritten each iteration). Useful for "
                             "identifying which cluster triggered a native crash in "
                             "pyfamsa/pytrimal/pyhmmer — Python try/except cannot catch "
                             "std::terminate() from C++ bindings, so if the process is "
                             "SIGABRT'd, this file holds the last attempted cluster id.")

    parser.add_argument("--trim", action="store_true",
                        help="Enable pytrimal alignment trimming before HMM construction. "
                             "Off by default. Trimming is not strictly necessary for HMM "
                             "building — hmmbuild's --symfrac already routes gappy columns "
                             "to insert states, and trimAl was benchmarked for ML tree "
                             "quality rather than HMM sensitivity. If your input clusters "
                             "are already well-separated by upstream clustering, trimming "
                             "usually adds little. It is also a known crash source for "
                             "certain malformed inputs (pytrimal can call std::terminate() "
                             "on unexpected residues, which Python cannot catch).")
    parser.add_argument("--trim_method", default="gappyout",
                        help="pytrimal AutomaticTrimmer method when --trim is set. "
                             "'gappyout' (default) uses only gap-distribution statistics "
                             "and never touches the similarity matrix, so it is robust to "
                             "residues the matrix doesn't recognize. 'automated1' also "
                             "incorporates similarity scores but can throw std::out_of_range "
                             "on certain inputs. Other methods from pytrimal are accepted "
                             "verbatim.")
    parser.add_argument("--min_seqs_for_trim", type=int, default=3,
                        help="Minimum sequences in cluster to apply trimming (default: 3)")

    parser.add_argument("--include_singletons", action="store_true",
                        help="Include singleton clusters (n=1) as single-sequence HMMs "
                             "(phmmer-equivalent), via pyhmmer's Builder.build(). "
                             "Default: exclude them. Single-sequence HMMs cannot capture "
                             "family-level conservation, their calibration is degenerate "
                             "(TC = self-score represents training-set memorization, not a "
                             "typical TP score), and they inflate the multiple-testing burden "
                             "with models that fire on essentially one sequence each. When "
                             "included, singletons are always calibrated with the 'simple' "
                             "method regardless of --calibration_method, since k-fold CV is "
                             "undefined for n=1.")

    parser.add_argument("--no_calibrate", action="store_true",
                        help="Disable cutoff calibration (default: on)")
    parser.add_argument("--calibration_method", choices=["kfold", "simple"], default="simple",
                        help="Calibration method for multi-member clusters. "
                             "'simple' = max-FP / min-TP from single hmmsearch (fast, default). "
                             "'kfold' = k-fold CV + F-measure (slow but more rigorous). "
                             "Singletons always use 'simple'.")
    parser.add_argument("--k_folds", type=int, default=3,
                        help="Number of CV folds for k-fold calibration (default: 3)")
    parser.add_argument("--epsilon", type=float, default=0.1,
                        help="Bit score margin above NC when setting GA in simple mode (default: 0.1)")
    parser.add_argument("--force_score_members", action="store_true",
                        help="Force-score cluster members that fall below hmmsearch inclusion "
                             "thresholds so they contribute to the TP distribution")

    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Threads for pyfamsa/pyhmmer (default: 1)")
    parser.add_argument("--random_state", type=int, default=42,
                        help="Random state for reproducible CV fold assignment (default: 42)")
    parser.add_argument("--log_level", default="INFO", help="loguru log level (default: INFO)")
    opts = parser.parse_args()

    logger.remove()
    logger.add(sys.stderr, level=opts.log_level)

    alphabet = Alphabet.amino()
    background = Background(alphabet)
    builder = Builder(alphabet, seed=opts.random_state)

    if opts.msa_directory is not None:
        Path(opts.msa_directory).mkdir(parents=True, exist_ok=True)

    protein_to_seq = load_proteins(opts.fasta)
    cluster_to_proteins = load_clusters(opts.clusters)

    cluster_to_hmm, singleton_cluster_ids = build_all_hmms(
        cluster_to_proteins=cluster_to_proteins,
        protein_to_seq=protein_to_seq,
        do_trim=opts.trim,
        trim_method=opts.trim_method,
        min_seqs_for_trim=opts.min_seqs_for_trim,
        alphabet=alphabet,
        builder=builder,
        background=background,
        threads=opts.threads,
        include_singletons=opts.include_singletons,
        msa_directory=opts.msa_directory,
        progress_file=opts.progress_file,
    )

    if not cluster_to_hmm:
        logger.error("No HMMs built, exiting")
        sys.exit(1)

    if not opts.no_calibrate:
        logger.info(f"Calibrating cutoffs using method={opts.calibration_method}")
        digital_proteome = digitize_proteome(protein_to_seq, alphabet)
        # For simple calibration we only need bit scores from hmmsearch; the
        # raw text sequences are no longer used. For kfold we must keep them
        # because CV rebuilds HMMs on each fold. Freeing the text proteome
        # matters for large inputs (e.g. ~6.8M proteins) where it can be
        # several GB on top of the digital proteome.
        if opts.calibration_method == "simple":
            n_freed = len(protein_to_seq)
            protein_to_seq_for_calibrate = None
            del protein_to_seq
            gc.collect()
            logger.info(f"Freed text proteome ({n_freed} sequences) prior to simple calibration")
        else:
            protein_to_seq_for_calibrate = protein_to_seq
        df_stats = calibrate_all(
            method=opts.calibration_method,
            cluster_to_hmm=cluster_to_hmm,
            cluster_to_proteins=cluster_to_proteins,
            protein_to_seq=protein_to_seq_for_calibrate,
            digital_proteome=digital_proteome,
            alphabet=alphabet,
            builder=builder,
            background=background,
            do_trim=opts.trim,
            trim_method=opts.trim_method,
            min_seqs_for_trim=opts.min_seqs_for_trim,
            k_folds=opts.k_folds,
            epsilon=opts.epsilon,
            force_score_members=opts.force_score_members,
            threads=opts.threads,
            random_state=opts.random_state,
            singleton_cluster_ids=singleton_cluster_ids,
        )
        if opts.output_stats:
            df_stats.to_csv(opts.output_stats, sep="\t", index=False)
            logger.info(f"Wrote cutoff stats to {opts.output_stats}")
        logger.info(f"Calibration summary:\n{df_stats['note'].value_counts().to_string()}")

    write_hmms(cluster_to_hmm, opts.output_hmm)
    dest = "stdout" if opts.output_hmm == "-" else opts.output_hmm
    logger.info(f"Wrote {len(cluster_to_hmm)} HMMs to {dest}")


if __name__ == "__main__":
    main()
