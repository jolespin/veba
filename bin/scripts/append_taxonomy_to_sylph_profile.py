#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, gzip, logging
import pandas as pd

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2025.1.20"

# Setup logger
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

def load_taxonomy_mapping(filepath):
    """
    Load taxonomy mapping from TSV or TSV.GZ file.
    Returns dict: {genome_key: taxonomy_string}
    """
    mapping = {}
    
    if filepath.endswith('.gz'):
        open_func = gzip.open
        mode = 'rt'
    else:
        open_func = open
        mode = 'r'
    
    with open_func(filepath, mode) as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    genome_key = parts[0]
                    taxonomy = parts[1]
                    mapping[genome_key] = taxonomy
    
    return mapping

def main(args=None):
    # Path info
    script_directory = os.path.dirname(os.path.abspath(__file__))
    script_filename = __program__
    
    description = """
    Running: {} v{} via Python v{} | {}
    
    Append taxonomy and database information to Sylph profiling results.
    Taxonomy mappings are matched against genome keys from the Genome_file column.
    If not found, falls back to Contig_name column.
    """.format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    
    usage = "{} -i <sylph_results.tsv> -t <tax1.tsv> <tax2.tsv> -n <db1> <db2> -o <output.tsv>".format(__program__)
    epilog = "Copyright 2025 Josh L. Espinoza (jespinoz@jcvi.org)"
    
    # Parser
    parser = argparse.ArgumentParser(
        description=description, 
        usage=usage, 
        epilog=epilog, 
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument(
        "-i", "--sylph_results",
        type=str,
        required=True,
        help="path/to/sylph_results.tsv"
    )
    parser.add_argument(
        "-t", "--taxonomy_mappings",
        type=str,
        nargs='+',
        required=True,
        help="One or more taxonomy mapping files (TSV or TSV.GZ, no header, genome_key<tab>taxonomy)"
    )
    parser.add_argument(
        "-n", "--database_names",
        type=str,
        nargs='+',
        required=True,
        help="Database names corresponding to taxonomy mappings (same order, same count)"
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default="stdout",
        help="Output filepath (.tsv, .tsv.gz, or .parquet) [Default: stdout]"
    )
    
    opts = parser.parse_args()
    opts.script_directory = script_directory
    opts.script_filename = script_filename
    
    # Validate inputs
    if len(opts.taxonomy_mappings) != len(opts.database_names):
        logger.error(
            f"Number of taxonomy mappings ({len(opts.taxonomy_mappings)}) "
            f"does not match number of database names ({len(opts.database_names)})"
        )
        sys.exit(1)
    
    # Check that taxonomy files exist
    for tax_file in opts.taxonomy_mappings:
        if not os.path.exists(tax_file):
            logger.error(f"Taxonomy mapping file does not exist: {tax_file}")
            sys.exit(1)
    
    # Load taxonomy mappings and create dictionaries
    logger.info(f"Loading {len(opts.taxonomy_mappings)} taxonomy mapping file(s)...")
    
    genome_to_taxonomy = {}
    genome_to_database = {}
    
    for db_name, tax_file in zip(opts.database_names, opts.taxonomy_mappings):
        logger.info(f"Loading {db_name}: {tax_file}")
        tax_mapping = load_taxonomy_mapping(tax_file)
        
        for genome_key, taxonomy in tax_mapping.items():
            genome_to_taxonomy[genome_key] = taxonomy
            genome_to_database[genome_key] = db_name
        
        logger.info(f"  Loaded {len(tax_mapping)} genome-taxonomy mappings")
    
    logger.info(f"Total unique genomes with taxonomy: {len(genome_to_taxonomy)}")
    
    # Load Sylph results
    logger.info(f"Loading Sylph results: {opts.sylph_results}")
    df = pd.read_csv(opts.sylph_results, sep="\t")
    logger.info(f"Loaded {df.shape[0]} rows, {df.shape[1]} columns")
    
    # Check required columns
    required_cols = ['Genome_file', 'Contig_name']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        logger.error(f"Missing required columns: {missing_cols}")
        logger.error(f"Available columns: {list(df.columns)}")
        sys.exit(1)
    
    # Add Taxonomy and Database columns
    logger.info("Mapping taxonomy and database information...")
    
    taxonomies = []
    databases = []
    not_found = []
    
    for idx, row in df.iterrows():
        genome_file = row['Genome_file']
        contig_name = row['Contig_name']
        
        # Try to find taxonomy using the Genome_file value first
        taxonomy = genome_to_taxonomy.get(genome_file)
        database = genome_to_database.get(genome_file)
        
        # If not found, try Contig_name
        if taxonomy is None:
            taxonomy = genome_to_taxonomy.get(contig_name)
            database = genome_to_database.get(contig_name)
        
        # Track if not found
        if taxonomy is None:
            not_found.append(genome_file)
        
        taxonomies.append(taxonomy)
        databases.append(database)
    
    # If any not found, fail immediately
    if not_found:
        unique_not_found = sorted(set(not_found))
        logger.error(f"FATAL: Could not find taxonomy for {len(not_found)} genome key(s)")
        logger.error(f"Unique genome keys without taxonomy: {len(unique_not_found)}")
        
        if len(unique_not_found) < 100:
            for key in unique_not_found:
                logger.error(f"  - {key}")
        
        logger.error("All genome keys must have taxonomy annotations. Exiting.")
        sys.exit(1)
    
    logger.info(f"Successfully mapped taxonomy for all {len(taxonomies)} entries")
    
    # Add new columns
    df['Database'] = databases
    df['Taxonomy'] = taxonomies
    
    # Write output
    if opts.output == "stdout":
        df.to_csv(sys.stdout, sep="\t", index=False)
    elif opts.output.endswith('.parquet'):
        logger.info(f"Writing output to: {opts.output} (parquet format)")
        df.to_parquet(opts.output, index=False)
    elif opts.output.endswith('.tsv.gz'):
        logger.info(f"Writing output to: {opts.output} (compressed TSV)")
        df.to_csv(opts.output, sep="\t", index=False, compression='gzip')
    elif opts.output.endswith('.tsv'):
        logger.info(f"Writing output to: {opts.output} (TSV)")
        df.to_csv(opts.output, sep="\t", index=False)
    else:
        logger.warning(f"Unknown output format for {opts.output}, writing as TSV")
        df.to_csv(opts.output, sep="\t", index=False)
    
    logger.info("Done!")

if __name__ == "__main__":
    main()