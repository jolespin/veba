#!/usr/bin/env python
"""
Build per-cluster HMMs from a protein fasta + cluster membership table.

Workspace model: writes everything under {output_directory}/, with
intermediate artifacts checkpointed so crashed runs can resume with the
same command.

Calibration methods:
  - 'simple' (default): NC = max bit score on non-members, TC = min bit
    score on members. GA placed between NC and TC. A single batched
    hmmsearch over all HMMs × the proteome.
  - 'kfold': k-fold CV with F-measure maximization (slow, per-cluster).
    Build HMM from k-1 folds, score held-out fold. GA = mean per-fold T.

Singletons are included by default as single-sequence HMMs (opt out with
--exclude_singletons). Singletons always use 'simple' calibration.
"""
import argparse
import gzip
import hashlib
import json
import os
import re
import shutil
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


# Valid fasta extensions for database name auto-detection
_FASTA_EXTENSIONS = (".fasta", ".fa", ".faa", ".fsa")


def infer_database_name(fasta_path):
    """Strip fasta extension (and optional .gz) from filename to get a database name.

    Returns None if the filename doesn't end with a recognized fasta extension.
    """
    name = os.path.basename(fasta_path)
    if name.endswith(".gz"):
        name = name[:-3]
    for ext in _FASTA_EXTENSIONS:
        if name.endswith(ext):
            return name[: -len(ext)]
    return None


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
        name=cluster_id,
        sequences=[
            TextSequence(name=pid, sequence=seq)
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
    text_seq = TextSequence(name=protein_id, sequence=seq_str)
    digital_seq = text_seq.digitize(alphabet)
    hmm, _, _ = builder.build(digital_seq, background)
    hmm.name = cluster_id
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
        ts = TextSequence(name=pid, sequence=seq)
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
    score_dict = {hit.name: hit.score for hit in top_hits}

    if force_members:
        missing = [m for m in force_members if m not in score_dict]
        if missing:
            missing_set = set(missing)
            subset = [ds for ds in digital_proteome if ds.name in missing_set]
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
                        name = hit.name
                        if name not in score_dict:
                            score_dict[name] = hit.score
                except Exception as e:
                    logger.debug(f"Force-scoring failed for HMM {hmm.name}: {e}")
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
    """
    h = hashlib.md5(f"{random_state}::{cluster_id}".encode()).hexdigest()
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
# Checkpointing
# ---------------------------------------------------------------------------

def checkpoint_uncalibrated_hmm_path(intermediate_directory, cluster_id):
    """Path where an uncalibrated HMM is stored (Phase 1 output)."""
    safe = sanitize_for_filename(cluster_id)
    return Path(intermediate_directory) / "uncalibrated_hmms" / f"{safe}.hmm.gz"


def checkpoint_calibrated_hmm_path(intermediate_directory, cluster_id):
    """Path where a calibrated HMM is stored (Phase 2 output)."""
    safe = sanitize_for_filename(cluster_id)
    return Path(intermediate_directory) / "calibrated_hmms" / f"{safe}.hmm.gz"


def checkpoint_stats_path(intermediate_directory, cluster_id):
    """Path where the stats row for a cluster is stored, one JSON object per file."""
    safe = sanitize_for_filename(cluster_id)
    return Path(intermediate_directory) / "stats" / f"{safe}.json"


def checkpoint_singleton_marker_path(intermediate_directory, cluster_id):
    """Empty file marking a cluster as singleton (preserved across resume)."""
    safe = sanitize_for_filename(cluster_id)
    return Path(intermediate_directory) / "singleton_markers" / f"{safe}.marker"


def write_hmm_gzip_atomic(hmm, path_output):
    """Write HMM to {path}.tmp then rename to {path}. Safe against mid-write crashes.

    A resume-time "does this cluster have a checkpoint" check looks for the
    final filename; the .tmp file is invisible to that check.
    """
    path_output = Path(path_output)
    path_tmp = path_output.with_suffix(path_output.suffix + ".tmp")
    with gzip.open(path_tmp, "wb") as file_handle:
        hmm.write(file_handle)
    os.replace(path_tmp, path_output)


def write_stats_atomic(stats_row, path_output):
    """Write stats row as a one-line JSON atomically."""
    path_output = Path(path_output)
    path_tmp = path_output.with_suffix(path_output.suffix + ".tmp")
    with open(path_tmp, "w") as file_handle:
        print(json.dumps(stats_row), file=file_handle)
    os.replace(path_tmp, path_output)


def read_checkpoint_hmm(path_hmm):
    """Read a single HMM back from a gzipped checkpoint file."""
    with gzip.open(path_hmm, "rb") as file_handle:
        with pyhmmer.plan7.HMMFile(file_handle) as hmm_file:
            return next(iter(hmm_file))


def read_checkpoint_stats(path_stats):
    """Read a stats row back from JSON checkpoint."""
    with open(path_stats) as file_handle:
        return json.loads(file_handle.read())


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------


def _worker_build_cluster(task):
    """Subprocess worker: build one cluster's HMM and return it as bytes.

    Each worker creates its own Alphabet/Background/Builder to avoid sharing
    pyhmmer C-side state across processes. With maxtasksperchild=1, the worker
    dies after every task, so Easel memory state cannot accumulate across
    clusters — this is the whole reason for subprocess mode.

    task is a dict with: cluster_id, present_pairs (list of (pid, seq) tuples),
    is_singleton, do_trim, trim_method, min_seqs_for_trim, random_state,
    msa_directory (str or None).

    Returns (cluster_id, hmm_bytes_or_None, is_singleton, error_message).
    """
    # These imports happen in the worker process. They're cheap because Python
    # caches them after the first task of this worker's lifetime, but with
    # maxtasksperchild=1 we pay the import cost per cluster. It's ~100ms total.
    from pyhmmer.easel import Alphabet as _Alphabet
    from pyhmmer.plan7 import Builder as _Builder, Background as _Background
    import io

    cluster_id = task["cluster_id"]
    present_pairs = task["present_pairs"]
    is_singleton = task["is_singleton"]
    msa_directory = task["msa_directory"]

    try:
        alphabet = _Alphabet.amino()
        background = _Background(alphabet)
        builder = _Builder(alphabet, seed=task["random_state"])

        if is_singleton:
            pid, seq_str = present_pairs[0]
            hmm = build_hmm_from_single_sequence(
                cluster_id, pid, seq_str, alphabet, builder, background,
            )
            if msa_directory is not None:
                try:
                    write_msa_gzip(cluster_id, [(pid, seq_str)], msa_directory)
                except Exception as e:
                    # Non-fatal
                    pass
        else:
            hmm = make_hmm_for_cluster(
                cluster_id=cluster_id,
                protein_ids=[pid for pid, _ in present_pairs],
                protein_to_seq={pid: seq for pid, seq in present_pairs},
                do_trim=task["do_trim"],
                trim_method=task["trim_method"],
                min_seqs_for_trim=task["min_seqs_for_trim"],
                alphabet=alphabet,
                builder=builder,
                background=background,
                threads=1,  # pyfamsa threads=1 per worker; parallelism is across workers
                msa_directory=msa_directory,
            )
            if hmm is None:
                return (cluster_id, None, is_singleton, "build returned None")

        # Serialize HMM to bytes for IPC
        buf = io.BytesIO()
        hmm.write(buf)
        return (cluster_id, buf.getvalue(), is_singleton, None)
    except Exception as e:
        return (cluster_id, None, is_singleton, f"{type(e).__name__}: {e}")


def phase1_build_all(
    cluster_to_proteins,
    protein_to_seq,
    do_trim,
    trim_method,
    min_seqs_for_trim,
    threads,
    exclude_singletons,
    intermediate_directory,
    msa_directory,
    random_state,
):
    """Phase 1: build uncalibrated HMMs for every cluster, in subprocesses.

    Each cluster is built in a fresh subprocess (maxtasksperchild=1) to prevent
    pyhmmer/Easel C-side memory state from accumulating across clusters. Up to
    `threads` clusters are built in parallel.

    HMMs are written to disk as bytes streamed back from workers. Resume: skip
    clusters whose uncalibrated HMM file already exists.

    Returns (singleton_cluster_ids, cluster_ids_with_hmms).
    """
    import multiprocessing
    from concurrent.futures import ProcessPoolExecutor, as_completed

    singleton_cluster_ids = set()
    cluster_ids_with_hmms = []

    # Build list of tasks, skipping already-checkpointed clusters
    tasks_to_submit = []
    n_resumed = 0
    for cluster_id, protein_ids in cluster_to_proteins.items():
        path_hmm = checkpoint_uncalibrated_hmm_path(intermediate_directory, cluster_id)
        if path_hmm.exists():
            cluster_ids_with_hmms.append(cluster_id)
            if checkpoint_singleton_marker_path(intermediate_directory, cluster_id).exists():
                singleton_cluster_ids.add(cluster_id)
            n_resumed += 1
            continue

        present = [p for p in protein_ids if p in protein_to_seq]
        if len(present) == 0:
            logger.warning(f"Cluster {cluster_id}: no members present in fasta — skipping")
            continue

        is_singleton = len(present) == 1
        if is_singleton and exclude_singletons:
            continue

        tasks_to_submit.append({
            "cluster_id": cluster_id,
            "present_pairs": [(p, protein_to_seq[p]) for p in present],
            "is_singleton": is_singleton,
            "do_trim": do_trim,
            "trim_method": trim_method,
            "min_seqs_for_trim": min_seqs_for_trim,
            "random_state": random_state,
            "msa_directory": str(msa_directory) if msa_directory is not None else None,
        })

    if n_resumed:
        logger.info(f"Phase 1: resumed {n_resumed} uncalibrated HMMs from checkpoint")

    n_built = 0
    n_skipped = 0

    if not tasks_to_submit:
        logger.info("Phase 1: nothing to build")
        return singleton_cluster_ids, cluster_ids_with_hmms

    logger.info(
        f"Phase 1: building {len(tasks_to_submit)} HMMs in subprocesses "
        f"({threads} parallel workers, fresh process per cluster)"
    )

    # ProcessPoolExecutor with maxtasksperchild=1 forces worker turnover after
    # every task. This is the entire point: Easel memory state is reset per
    # cluster, preventing accumulation that was causing SIGABRT in the main run.
    mp_context = multiprocessing.get_context("spawn")  # spawn is cleaner than fork on macOS
    with ProcessPoolExecutor(
        max_workers=threads,
        mp_context=mp_context,
        max_tasks_per_child=1,
    ) as executor:
        futures = {
            executor.submit(_worker_build_cluster, task): task["cluster_id"]
            for task in tasks_to_submit
        }

        for future in tqdm(
            as_completed(futures),
            desc="Phase 1: building HMMs",
            total=len(futures),
            unit="cluster",
            mininterval=1.0,
        ):
            cluster_id = futures[future]
            try:
                _cid, hmm_bytes, is_singleton, err = future.result()
            except Exception as e:
                logger.warning(f"Cluster {cluster_id}: worker crashed ({type(e).__name__}: {e})")
                n_skipped += 1
                continue

            if hmm_bytes is None:
                logger.warning(f"Cluster {cluster_id}: build failed ({err})")
                n_skipped += 1
                continue

            if is_singleton:
                singleton_cluster_ids.add(cluster_id)

            path_hmm = checkpoint_uncalibrated_hmm_path(intermediate_directory, cluster_id)
            try:
                # Write HMM bytes atomically (.tmp → rename)
                path_tmp = path_hmm.with_suffix(path_hmm.suffix + ".tmp")
                with gzip.open(path_tmp, "wb") as fh:
                    fh.write(hmm_bytes)
                os.replace(path_tmp, path_hmm)
                if is_singleton:
                    checkpoint_singleton_marker_path(intermediate_directory, cluster_id).touch()
                cluster_ids_with_hmms.append(cluster_id)
                n_built += 1
            except Exception as e:
                logger.warning(f"Cluster {cluster_id}: Phase 1 checkpoint write failed ({e})")

    logger.info(
        f"Phase 1 complete. Built: {n_built}, resumed: {n_resumed}, skipped: {n_skipped}"
    )
    return singleton_cluster_ids, cluster_ids_with_hmms


def phase2_calibrate_batched_simple(
    cluster_ids_all,
    cluster_to_proteins,
    digital_proteome,
    threads,
    epsilon,
    singleton_cluster_ids,
    intermediate_directory,
):
    """Phase 2 (simple mode): one batched hmmsearch against the full proteome.

    HMMs are read lazily from disk, calibrated, and written back. Already-
    calibrated clusters are skipped so resumed runs don't rescore them.

    Returns df_stats.
    """
    cluster_to_members = {c: set(p) for c, p in cluster_to_proteins.items()}

    stats_rows = []
    clusters_to_do = []
    n_resumed = 0

    for cluster_id in cluster_ids_all:
        path_calibrated = checkpoint_calibrated_hmm_path(intermediate_directory, cluster_id)
        path_stats = checkpoint_stats_path(intermediate_directory, cluster_id)
        if path_calibrated.exists() and path_stats.exists():
            try:
                stats_rows.append(read_checkpoint_stats(path_stats))
                n_resumed += 1
                continue
            except Exception as e:
                logger.warning(
                    f"Cluster {cluster_id}: calibrated checkpoint unreadable ({e}), rebuilding"
                )
        clusters_to_do.append(cluster_id)

    if n_resumed:
        logger.info(f"Phase 2: resumed {n_resumed} calibrated HMMs from checkpoint")

    if not clusters_to_do:
        logger.info("Phase 2: nothing to calibrate")
        return pd.DataFrame(stats_rows) if stats_rows else pd.DataFrame()

    # Generator: lazily read each uncalibrated HMM as pyhmmer asks for it
    def hmm_iter():
        for cluster_id in clusters_to_do:
            yield read_checkpoint_hmm(
                checkpoint_uncalibrated_hmm_path(intermediate_directory, cluster_id)
            )

    logger.info(
        f"Phase 2: batched hmmsearch over {len(clusters_to_do)} HMMs × "
        f"{len(digital_proteome)} proteins"
    )
    tophits_stream = pyhmmer.hmmsearch(hmm_iter(), digital_proteome, cpus=threads)

    for cluster_id, top_hits in tqdm(
        zip(clusters_to_do, tophits_stream),
        desc="Phase 2: batched hmmsearch + cutoffs",
        total=len(clusters_to_do),
        unit="hmm",
        mininterval=1.0,
    ):
        # Reload HMM to write cutoffs
        hmm = read_checkpoint_hmm(
            checkpoint_uncalibrated_hmm_path(intermediate_directory, cluster_id)
        )

        members = cluster_to_members[cluster_id]
        is_singleton = cluster_id in singleton_cluster_ids

        score_dict = {hit.name: hit.score for hit in top_hits}
        tp_scores = [s for name, s in score_dict.items() if name in members]
        fp_scores = [s for name, s in score_dict.items() if name not in members]

        nc = max(fp_scores) if fp_scores else None
        tc = min(tp_scores) if tp_scores else None

        if nc is None and tc is None:
            ga, tc_out, nc_out, note = 0.0, 0.0, 0.0, "no_hits"
        elif nc is None:
            ga = max(tc - epsilon, 0.0)
            tc_out, nc_out, note = tc, 0.0, "clean"
        elif tc is None:
            ga, tc_out, nc_out, note = nc + epsilon, nc + epsilon, nc, "no_self_hits"
        elif tc <= nc:
            ga, tc_out, nc_out, note = nc + epsilon, tc, nc, "overlap"
        else:
            ga = min(nc + epsilon, tc)
            tc_out, nc_out, note = tc, nc, "ok"

        if is_singleton:
            note = f"singleton_{note}"

        hmm.cutoffs.gathering = (ga, ga)
        hmm.cutoffs.trusted = (tc_out, tc_out)
        hmm.cutoffs.noise = (nc_out, nc_out)

        stats_row = {
            "id_protein_cluster": cluster_id,
            "n_members": len(members),
            "ga": ga, "tc": tc_out, "nc": nc_out,
            "n_tp_hits": len(tp_scores),
            "n_fp_hits": len(fp_scores),
            "note": note,
        }
        stats_rows.append(stats_row)

        try:
            write_hmm_gzip_atomic(
                hmm, checkpoint_calibrated_hmm_path(intermediate_directory, cluster_id)
            )
            write_stats_atomic(
                stats_row, checkpoint_stats_path(intermediate_directory, cluster_id)
            )
        except Exception as e:
            logger.warning(f"Cluster {cluster_id}: Phase 2 checkpoint write failed ({e})")
        del hmm

    return pd.DataFrame(stats_rows) if stats_rows else pd.DataFrame()


def phase2_calibrate_per_cluster_kfold(
    cluster_ids_all,
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
    threads,
    force_score_members,
    random_state,
    singleton_cluster_ids,
    intermediate_directory,
    epsilon,
):
    """Phase 2 (kfold mode): per-cluster k-fold CV calibration.

    Singletons fall back to simple calibration (kfold undefined for n=1).

    Returns df_stats.
    """
    cluster_to_members = {c: set(p) for c, p in cluster_to_proteins.items()}
    stats_rows = []
    n_resumed = 0

    clusters_to_do = []
    for cluster_id in cluster_ids_all:
        path_calibrated = checkpoint_calibrated_hmm_path(intermediate_directory, cluster_id)
        path_stats = checkpoint_stats_path(intermediate_directory, cluster_id)
        if path_calibrated.exists() and path_stats.exists():
            try:
                stats_rows.append(read_checkpoint_stats(path_stats))
                n_resumed += 1
                continue
            except Exception as e:
                logger.warning(
                    f"Cluster {cluster_id}: calibrated checkpoint unreadable ({e}), rebuilding"
                )
        clusters_to_do.append(cluster_id)

    if n_resumed:
        logger.info(f"Phase 2: resumed {n_resumed} calibrated HMMs from checkpoint")

    for cluster_id in tqdm(
        clusters_to_do,
        desc="Phase 2: kfold calibration",
        unit="cluster",
        mininterval=1.0,
    ):
        hmm = read_checkpoint_hmm(
            checkpoint_uncalibrated_hmm_path(intermediate_directory, cluster_id)
        )

        members = cluster_to_members[cluster_id]
        is_singleton = cluster_id in singleton_cluster_ids

        if is_singleton:
            result = calibrate_simple(
                cluster_id=cluster_id,
                hmm=hmm,
                members_set=members,
                digital_proteome=digital_proteome,
                threads=threads,
                epsilon=epsilon,
                force_score_members=False,
                alphabet=alphabet,
                background=background,
            )
            result["note"] = f"singleton_{result['note']}"
        else:
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

        ga, tc, nc = result["ga"], result["tc"], result["nc"]
        hmm.cutoffs.gathering = (ga, ga)
        hmm.cutoffs.trusted = (tc, tc)
        hmm.cutoffs.noise = (nc, nc)

        stats_row = {
            "id_protein_cluster": cluster_id,
            "n_members": len(members),
            "ga": ga, "tc": tc, "nc": nc,
            "n_tp_hits": result["n_tp"],
            "n_fp_hits": result["n_fp"],
            "note": result["note"],
        }
        if not is_singleton:
            stats_row["cv_thresholds"] = ",".join(
                f"{t:.3f}" for t in result.get("cv_thresholds", [])
            )
            stats_row["cv_fmeasure"] = ",".join(
                f"{f:.3f}" for f in result.get("cv_fmeasure", [])
            )
        stats_rows.append(stats_row)

        try:
            write_hmm_gzip_atomic(
                hmm, checkpoint_calibrated_hmm_path(intermediate_directory, cluster_id)
            )
            write_stats_atomic(
                stats_row, checkpoint_stats_path(intermediate_directory, cluster_id)
            )
        except Exception as e:
            logger.warning(f"Cluster {cluster_id}: Phase 2 checkpoint write failed ({e})")
        del hmm

    return pd.DataFrame(stats_rows) if stats_rows else pd.DataFrame()


def write_hmm_database(cluster_ids, intermediate_directory, output_path, calibrated):
    """Concatenate per-cluster HMMs from disk into a single gzipped HMM database.

    Streams one HMM at a time; never holds more than one HMM in memory.
    """
    subdir = "calibrated_hmms" if calibrated else "uncalibrated_hmms"
    path_intermediate = Path(intermediate_directory)

    path_output = Path(output_path)
    path_tmp = path_output.with_suffix(path_output.suffix + ".tmp")

    with gzip.open(path_tmp, "wb") as out_handle:
        for cluster_id in tqdm(
            cluster_ids, desc="Writing HMM database", unit="hmm", mininterval=1.0
        ):
            path_hmm = path_intermediate / subdir / f"{sanitize_for_filename(cluster_id)}.hmm.gz"
            if not path_hmm.exists():
                # Fallback: if calibration was run but this cluster didn't finish, use uncalibrated
                if calibrated:
                    fallback = path_intermediate / "uncalibrated_hmms" / f"{sanitize_for_filename(cluster_id)}.hmm.gz"
                    if fallback.exists():
                        path_hmm = fallback
                    else:
                        logger.warning(f"Cluster {cluster_id}: HMM missing on disk, skipping")
                        continue
                else:
                    logger.warning(f"Cluster {cluster_id}: HMM missing on disk, skipping")
                    continue

            hmm = read_checkpoint_hmm(path_hmm)
            hmm.write(out_handle)
            del hmm

    os.replace(path_tmp, path_output)

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build per-cluster HMMs from a protein fasta and cluster table. "
            "Everything is written under an output directory, with intermediate "
            "artifacts checkpointed so interrupted runs can resume with the same command."
        ),
    )
    parser.add_argument("-f", "--fasta", required=True, help="Protein fasta file")
    parser.add_argument("-c", "--clusters", required=True,
                        help="TSV (no header): id_protein\\tid_protein_cluster")
    parser.add_argument("-o", "--output_directory", required=True,
                        help="Output directory. Will contain the final HMM database, "
                             "stats TSV, and an intermediate/ subdirectory with "
                             "per-cluster checkpoints.")
    parser.add_argument("-n", "--database_name", default=None,
                        help="Database name, used as prefix for output files "
                             "({output_directory}/{database_name}.hmm.gz, .stats.tsv). "
                             "If not specified, auto-detected from the fasta filename "
                             "(strips .fasta/.fa/.faa/.fsa and optional .gz).")
    parser.add_argument("--force", action="store_true",
                        help="Wipe the intermediate directory before running, forcing a "
                             "full rebuild. Default: resume from existing checkpoints.")
    parser.add_argument("--remove_intermediate_msa", action="store_true",
                        help="After a successful run, delete the intermediate MSA "
                             "directory to reclaim disk. Default: keep MSAs.")

    parser.add_argument("--no_trim", action="store_true",
                        help="Disable pytrimal trimming (default: on)")
    parser.add_argument("--trim_method", default="automated1",
                        help="pytrimal AutomaticTrimmer method (default: automated1)")
    parser.add_argument("--min_seqs_for_trim", type=int, default=3,
                        help="Minimum sequences in cluster to apply trimming (default: 3)")

    parser.add_argument("--exclude_singletons", action="store_true",
                        help="Exclude singleton clusters (n=1). Default: include them "
                             "as single-sequence HMMs (phmmer-equivalent).")

    parser.add_argument("--no_calibrate", action="store_true",
                        help="Disable cutoff calibration (default: on)")
    parser.add_argument("--calibration_method", choices=["kfold", "simple"], default="simple",
                        help="Calibration method for multi-member clusters. "
                             "'simple' = max-FP / min-TP from single batched hmmsearch (fast, default). "
                             "'kfold' = k-fold CV + F-measure (slow). "
                             "Singletons always use 'simple'.")
    parser.add_argument("--k_folds", type=int, default=3,
                        help="Number of CV folds for k-fold calibration (default: 3)")
    parser.add_argument("--epsilon", type=float, default=0.1,
                        help="Bit score margin above NC when setting GA (default: 0.1)")
    parser.add_argument("--force_score_members", action="store_true",
                        help="Force-score cluster members below hmmsearch inclusion "
                             "thresholds so they contribute to the TP distribution (kfold only)")

    parser.add_argument("--n_concurrent_hmm_workers", type=int, default=1,
                        help="Number of clusters to build in parallel during Phase 1. "
                             "Each cluster runs in a fresh subprocess to prevent Easel "
                             "memory state accumulation. (default: 1)")
    parser.add_argument("--n_threads_hmmsearch", type=int, default=1,
                        help="Threads for pyhmmer.hmmsearch during Phase 2 calibration "
                             "(default: 1)")
    parser.add_argument("--random_state", type=int, default=42,
                        help="Random state for reproducible CV fold assignment (default: 42)")
    parser.add_argument("--log_level", default="INFO", help="loguru log level (default: INFO)")
    opts = parser.parse_args()

    logger.remove()
    logger.add(sys.stderr, level=opts.log_level)

    # Resolve database name
    database_name = opts.database_name
    if database_name is None:
        database_name = infer_database_name(opts.fasta)
        if database_name is None:
            parser.error(
                f"Could not infer database name from fasta filename '{opts.fasta}'. "
                f"Expected extension {_FASTA_EXTENSIONS} (optionally .gz). "
                f"Please specify -n/--database_name explicitly."
            )
        logger.info(f"Auto-detected database name: {database_name}")

    # Set up workspace
    path_output = Path(opts.output_directory)
    path_intermediate = path_output / "intermediate"

    if opts.force and path_intermediate.exists():
        logger.info(f"--force: removing existing {path_intermediate}")
        shutil.rmtree(path_intermediate)

    (path_intermediate / "uncalibrated_hmms").mkdir(parents=True, exist_ok=True)
    (path_intermediate / "calibrated_hmms").mkdir(parents=True, exist_ok=True)
    (path_intermediate / "stats").mkdir(parents=True, exist_ok=True)
    (path_intermediate / "singleton_markers").mkdir(parents=True, exist_ok=True)
    msa_directory = path_intermediate / "msa"
    msa_directory.mkdir(parents=True, exist_ok=True)

    path_final_hmm = path_output / f"{database_name}.hmm.gz"
    path_final_stats = path_output / f"{database_name}.stats.tsv"

    # Core setup
    alphabet = Alphabet.amino()
    background = Background(alphabet)
    builder = Builder(alphabet, seed=opts.random_state)

    protein_to_seq = load_proteins(opts.fasta)
    cluster_to_proteins = load_clusters(opts.clusters)

    # ---- Phase 1: build uncalibrated HMMs (subprocess per cluster) ----
    singleton_cluster_ids, all_cluster_ids = phase1_build_all(
        cluster_to_proteins=cluster_to_proteins,
        protein_to_seq=protein_to_seq,
        do_trim=not opts.no_trim,
        trim_method=opts.trim_method,
        min_seqs_for_trim=opts.min_seqs_for_trim,
        threads=opts.n_concurrent_hmm_workers,
        exclude_singletons=opts.exclude_singletons,
        intermediate_directory=path_intermediate,
        msa_directory=msa_directory,
        random_state=opts.random_state,
    )

    if not all_cluster_ids:
        logger.error("No HMMs produced in Phase 1, exiting")
        sys.exit(1)

    # ---- Phase 2: calibrate (if enabled) ----
    df_stats = pd.DataFrame()
    if not opts.no_calibrate:
        logger.info("Digitizing proteome for calibration")
        digital_proteome = digitize_proteome(protein_to_seq, alphabet)

        if opts.calibration_method == "simple":
            del protein_to_seq  # not needed for simple
            df_stats = phase2_calibrate_batched_simple(
                cluster_ids_all=all_cluster_ids,
                cluster_to_proteins=cluster_to_proteins,
                digital_proteome=digital_proteome,
                threads=opts.n_threads_hmmsearch,
                epsilon=opts.epsilon,
                singleton_cluster_ids=singleton_cluster_ids,
                intermediate_directory=path_intermediate,
            )
        elif opts.calibration_method == "kfold":
            df_stats = phase2_calibrate_per_cluster_kfold(
                cluster_ids_all=all_cluster_ids,
                cluster_to_proteins=cluster_to_proteins,
                protein_to_seq=protein_to_seq,
                digital_proteome=digital_proteome,
                alphabet=alphabet,
                builder=builder,
                background=background,
                do_trim=not opts.no_trim,
                trim_method=opts.trim_method,
                min_seqs_for_trim=opts.min_seqs_for_trim,
                k_folds=opts.k_folds,
                threads=opts.n_threads_hmmsearch,
                force_score_members=opts.force_score_members,
                random_state=opts.random_state,
                singleton_cluster_ids=singleton_cluster_ids,
                intermediate_directory=path_intermediate,
                epsilon=opts.epsilon,
            )

        del digital_proteome

    # Write final outputs
    if not df_stats.empty:
        df_stats.to_csv(path_final_stats, sep="\t", index=False)
        logger.info(f"Wrote cutoff stats to {path_final_stats}")
    if not df_stats.empty and "note" in df_stats.columns:
        logger.info(f"Calibration summary:\n{df_stats['note'].value_counts().to_string()}")

    write_hmm_database(
        cluster_ids=all_cluster_ids,
        intermediate_directory=path_intermediate,
        output_path=path_final_hmm,
        calibrated=not opts.no_calibrate,
    )
    logger.info(f"Wrote {len(all_cluster_ids)} HMMs to {path_final_hmm}")

    if opts.remove_intermediate_msa:
        logger.info(f"--remove_intermediate_msa: removing {msa_directory}")
        shutil.rmtree(msa_directory, ignore_errors=True)


if __name__ == "__main__":
    main()
