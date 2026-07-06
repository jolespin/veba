#!/usr/bin/env python
import sys
import os
import gzip
import argparse
import logging
from pathlib import Path
from typing import TextIO, Optional, List, Dict, Tuple

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2026.1.2"

# Setup logging to stderr
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    stream=sys.stderr
)
logger = logging.getLogger(__name__)


def open_file(filepath: Path, mode: str = 'rt'):
    """Open file, handling both plain text and gzip compressed files."""
    if filepath.suffix == '.gz':
        return gzip.open(filepath, mode)
    return open(filepath, mode)


def extract_cluster_name(filepath: Path) -> str:
    """Extract cluster name from filename (remove .txt or .txt.gz)."""
    name = filepath.name
    if name.endswith('.txt.gz'):
        return name[:-7]
    elif name.endswith('.txt'):
        return name[:-4]
    return name


def parse_hit_header(lines: List[str], start_idx: int) -> Optional[Tuple[Dict[str, str], int]]:
    """
    Parse a single hit header section.
    Returns (header_dict, next_line_index) or None if no valid hit found.
    """
    header = {}
    idx = start_idx
    
    # Skip empty lines
    while idx < len(lines) and not lines[idx].strip():
        idx += 1
    
    if idx >= len(lines):
        return None
    
    # Check if this is a numbered hit (e.g., "1. NZ_CP119108")
    line = lines[idx].strip()
    if not line or not line[0].isdigit():
        return None
    
    # Parse accession (first line after number)
    parts = line.split(None, 1)
    if len(parts) < 2:
        return None
    
    header['accession'] = parts[1].strip()
    idx += 1
    
    # Parse Source, Type, Number of proteins, Cumulative BLAST score
    while idx < len(lines):
        line = lines[idx].strip()
        
        if not line:
            idx += 1
            continue
        
        if line.startswith('Source:'):
            header['source'] = line.split('Source:', 1)[1].strip()
        elif line.startswith('Type:'):
            header['type'] = line.split('Type:', 1)[1].strip()
        elif line.startswith('Number of proteins with BLAST hits to this cluster:'):
            header['number_of_hits'] = line.split(':', 1)[1].strip()
        elif line.startswith('Cumulative BLAST score:'):
            header['cumulative_blast_score'] = line.split(':', 1)[1].strip()
        elif line.startswith('Table of genes'):
            # Found start of genes table, skip it
            idx += 1
            break
        elif line.startswith('Table of Blast hits'):
            # Found blast hits table header
            idx += 1
            break
        
        idx += 1
    
    return header, idx


def parse_blast_hits(lines: List[str], start_idx: int) -> Tuple[List[Dict[str, str]], int]:
    """
    Parse blast hits table starting from start_idx.
    Returns (list of hit dicts, next_line_index).
    """
    hits = []
    idx = start_idx
    
    # Skip until we find the blast hits table or hit end markers
    while idx < len(lines):
        line = lines[idx].strip()
        
        # Stop conditions
        if not line:
            idx += 1
            continue
        
        # Check if we've hit the next section marker
        if line.startswith('>>') or (line and line[0].isdigit() and '. ' in line[:5]):
            break
        
        # Skip header lines
        if line.startswith('Table of') or line.startswith('query gene'):
            idx += 1
            continue
        
        # Parse blast hit line
        # Format: query_gene\tsubject_gene\t%identity\tblast_score\t%coverage\te-value
        parts = line.split('\t')
        if len(parts) == 6:
            hit = {
                'query_gene': parts[0].strip(),
                'subject_gene': parts[1].strip(),
                'percent_identity': parts[2].strip(),
                'blast_score': parts[3].strip(),
                'percent_coverage': parts[4].strip(),
                'e_value': parts[5].strip()
            }
            hits.append(hit)
        
        idx += 1
    
    return hits, idx


def parse_clusterblast_file(filepath: Path, cluster_method: str) -> List[Dict[str, str]]:
    """
    Parse a single clusterblast file and return list of row dictionaries.
    """
    cluster_name = extract_cluster_name(filepath)
    
    try:
        with open_file(filepath) as f:
            content = f.read()
    except Exception as e:
        logger.critical(f"Failed to read file {filepath}: {e}")
        return []
    
    if not content.strip():
        logger.critical(f"Empty file encountered: {filepath}")
        return []
    
    lines = content.split('\n')
    
    # Find "Significant hits:" section
    sig_hits_idx = None
    for i, line in enumerate(lines):
        if line.strip().startswith('Significant hits:'):
            sig_hits_idx = i
            break
    
    if sig_hits_idx is None:
        logger.warning(f"No 'Significant hits:' section found in {filepath}")
        return []
    
    # Check if there are any hits after "Significant hits:"
    has_hits = False
    for i in range(sig_hits_idx + 1, len(lines)):
        line = lines[i].strip()
        if line and line[0].isdigit():
            has_hits = True
            break
        if line.startswith('Details:'):
            break
    
    if not has_hits:
        logger.warning(f"No hits found in {filepath}")
        return []
    
    # Find "Details:" section
    details_idx = None
    for i, line in enumerate(lines):
        if line.strip().startswith('Details:'):
            details_idx = i
            break
    
    if details_idx is None:
        logger.warning(f"No 'Details:' section found in {filepath}")
        return []
    
    # Parse all hits in Details section
    all_rows = []
    idx = details_idx + 1
    
    while idx < len(lines):
        # Skip to next ">>" marker
        while idx < len(lines) and not lines[idx].strip().startswith('>>'):
            idx += 1
        
        if idx >= len(lines):
            break
        
        idx += 1  # Skip the ">>" line
        
        # Parse hit header
        result = parse_hit_header(lines, idx)
        if result is None:
            break
        
        header, idx = result
        
        # Parse blast hits for this header
        blast_hits, idx = parse_blast_hits(lines, idx)
        
        # Combine header with each blast hit
        for hit in blast_hits:
            row = {
                'cluster_method': cluster_method,
                'cluster': cluster_name,
                'accession': header.get('accession', ''),
                'source': header.get('source', ''),
                'type': header.get('type', ''),
                'number_of_hits': header.get('number_of_hits', ''),
                'cumulative_blast_score': header.get('cumulative_blast_score', ''),
                'query_gene': hit['query_gene'],
                'subject_gene': hit['subject_gene'],
                'percent_identity': hit['percent_identity'],
                'blast_score': hit['blast_score'],
                'percent_coverage': hit['percent_coverage'],
                'e_value': hit['e_value']
            }
            all_rows.append(row)
    
    return all_rows


def process_antismash_results(antismash_dir: Path) -> List[Dict[str, str]]:
    """
    Process all clusterblast subdirectories in antiSMASH results.
    """
    all_rows = []
    
    subdirs = ['clusterblast', 'knownclusterblast', 'subclusterblast']
    
    for subdir in subdirs:
        subdir_path = antismash_dir / subdir
        
        if not subdir_path.exists():
            logger.info(f"Subdirectory {subdir} not found, skipping")
            continue
        
        if not subdir_path.is_dir():
            logger.warning(f"{subdir_path} is not a directory, skipping")
            continue
        
        logger.info(f"Processing {subdir} directory")
        
        # Find all .txt and .txt.gz files
        txt_files = list(subdir_path.glob('*.txt'))
        txt_gz_files = list(subdir_path.glob('*.txt.gz'))
        all_files = txt_files + txt_gz_files
        
        logger.info(f"Found {len(all_files)} files in {subdir}")
        
        for filepath in all_files:
            logger.info(f"Parsing {filepath.name}")
            rows = parse_clusterblast_file(filepath, subdir)
            all_rows.extend(rows)
            logger.info(f"Extracted {len(rows)} rows from {filepath.name}")
    
    return all_rows


def write_output(rows: List[Dict[str, str]], output_path: Optional[str]):
    """
    Write rows to output file or stdout.
    """
    columns = [
        'cluster_method',
        'cluster',
        'accession',
        'source',
        'type',
        'number_of_hits',
        'cumulative_blast_score',
        'query_gene',
        'subject_gene',
        'percent_identity',
        'blast_score',
        'percent_coverage',
        'e_value'
    ]
    
    # Determine output handle
    if output_path is None or output_path == '-' or output_path == 'stdout':
        out_handle = sys.stdout
    else:
        output_path_obj = Path(output_path)
        if output_path_obj.suffix == '.gz':
            out_handle = gzip.open(output_path, 'wt')
        else:
            out_handle = open(output_path, 'w')
    
    try:
        # Write header
        print('\t'.join(columns), file=out_handle)
        
        # Write rows
        for row in rows:
            values = [row.get(col, '') for col in columns]
            print('\t'.join(values), file=out_handle)
        
        logger.info(f"Wrote {len(rows)} total rows to output")
    
    finally:
        if out_handle not in (sys.stdout, sys.stderr):
            out_handle.close()


def main():
    usage = "{} -i <antismash_results_directory> -o <output.tsv[.gz]>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"
    
    parser = argparse.ArgumentParser(
        description='Parse antiSMASH clusterblast results',
        usage=usage,
        epilog=epilog,
    )
    parser.add_argument(
        '-i', '--antismash_results',
        required=True,
        type=str,
        help='Path to antiSMASH results directory'
    )
    parser.add_argument(
        '-o', '--output',
        type=str,
        default=None,
        help='Output TSV file (default: stdout). Can be .gz compressed'
    )
    
    opts = parser.parse_args()
    
    antismash_dir = Path(opts.antismash_results)
    
    if not antismash_dir.exists():
        logger.error(f"antiSMASH results directory does not exist: {antismash_dir}")
        sys.exit(1)
    
    if not antismash_dir.is_dir():
        logger.error(f"Path is not a directory: {antismash_dir}")
        sys.exit(1)
    
    logger.info(f"Processing antiSMASH results from: {antismash_dir}")
    
    rows = process_antismash_results(antismash_dir)
    
    write_output(rows, opts.output)
    
    logger.info("Processing complete")


if __name__ == '__main__':
    main()