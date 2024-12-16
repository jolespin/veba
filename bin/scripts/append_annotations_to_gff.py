#!/usr/bin/env python
import sys, os, glob, argparse, gzip
from collections import OrderedDict, defaultdict
import pandas as pd
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2024.11.11"

def get_annotations_from_veba(f):
    """
    Parse a tab-delimited annotation file and extract annotations for each protein.

    Parameters
    ----------
    f : file-like object
        A file handle to the annotation file, expected to be in a specific tab-delimited format.

    Returns
    -------
    defaultdict
        A dictionary where each key is a protein identifier and each value is an OrderedDict 
        containing the extracted annotations per database including product names.

    Notes
    -----
    This function processes annotations from multiple databases such as UniRef, MIBiG, VFDB, 
    and CAZy for sseqid and Pfam, KOfam, NCBIfam-AMR, AntiFam for ids.
    """
    id_to_annotation = defaultdict(OrderedDict)

    df_annotations = pd.read_csv(f, sep="\t", index_col=0, header=[0,1])
    for id, row in tqdm(df_annotations.iterrows(), desc="Parsing annotations", unit=" annotations"):
        id_to_annotation[id]["product"] = row[("Consensus", "composite_name")]
        
        for database_name in ["UniRef", "MIBiG", "VFDB", "CAZy"]:
            annotation = row[(database_name, "sseqid")]
            if pd.notnull(annotation):
                id_to_annotation[id][database_name] = annotation
        for database_name in ["Pfam","KOfam", "NCBIfam-AMR", "AntiFam", "Enzymes"]:
            annotation = eval(row[(database_name, "ids")])
            if annotation:
                id_to_annotation[id][database_name] = ",".join(annotation)  
    return id_to_annotation

def get_annotations_from_pyhmmsearch(f, attribute_label):
    """
    Parse PyKOfamSearch/PyHMMerSearch output file to extract annotations for each protein
    
    Parameters
    ----------
    f : filehandle
        filehandle to the PyHMMerSearch output file
    attribute_label : str
        the label to use for the annotation attribute
    
    Returns
    -------
    id_to_annotation : dict
        dict of protein ids to annotation dictionaries
    """
    id_to_annotation = defaultdict(set)
    for line in tqdm(f, desc="Parsing annotations", unit=" annotations"):
        if not line.startswith("#"):
            id, annotation, *_ = line.split("\t")
            id_to_annotation[id].add(annotation)
    id_to_annotation = {k:{attribute_label:",".join(v)} for k,v in id_to_annotation.items()}
    return id_to_annotation
            
def get_annotations_from_pyhmmsearch_reformatted(f, attribute_label):
    
    """
    Parse a PyKOfamSearch/PyHMMerSearch reformatted output file to extract annotations for each protein

    Parameters
    ----------
    f : filehandle
        filehandle to the PyHMMerSearch reformatted output file
    attribute_label : str
        the label to use for the annotation attribute

    Returns
    -------
    id_to_annotation : dict
        dict of protein ids to annotation dictionaries
    """
    id_to_annotation = defaultdict(dict)
    for line in tqdm(f, desc="Parsing annotations", unit=" annotations"):
        if not line.startswith("#"):
            id, num, annotation, *_ = line.split("\t")
            annotation = eval(annotation)
            if annotation:
                id_to_annotation[id][attribute_label] = ",".join(annotation)
    return id_to_annotation

    

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i <input.gff> -a gene_id -o <output.gff>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--input", type=str, default="stdin", help = "path/to/gene_models.gff[.gz] [Default: stdin]")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/gene_models.updated.gff[.gz] [Default: stdout]")
    parser.add_argument("-a","--annotations", type=str, required=True, help = "path/to/annotations.tsv[.gz]")
    parser.add_argument("-g","--identifier_attribute", type=str, default="gene_id",  help = "Identifier attribute [Default: gene_id]")
    parser.add_argument("-f","--annotation_format", type=str, default="veba", 
                        choices={
                            # VEBA
                            "veba",
                            # PyKoFamSearch
                            "pykofamsearch", 
                            "pykofamsearch-no-header", 
                            "pykofamsearch-reformatted", 
                            "pykofamsearch-reformatted-no-header", 
                            # PyHMMSearch
                            "pyhmmsearch", 
                            "pyhmmsearch-no-header", 
                            "pyhmmsearch-reformatted",
                            "pyhmmsearch-reformatted-no-header",
                            }, 
                        help = "Annotation format. [Default: veba]")
    parser.add_argument("--features", type=str, default="gene,CDS,mRNA",  help = "Features to add annotations separated by commas [Default: gene,CDS,mRNA]")
    parser.add_argument("-l","--pyhmmsearch_attribute_label", type=str, default="hmm_id",  help = "Attribute label for PyHMMSearch (e.g., Pfam, TIGRFAM).  If PyKOfamSearch then KOfam is used. [Default: hmm_id]")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Input
    if opts.input == "stdin":
        f_in = sys.stdin 
    else:
        if opts.input.endswith(".gz"):
            f_in = gzip.open(opts.input, "rt")
        else:
            f_in = open(opts.input, "r")

    # Output
    if opts.output == "stdout":
        f_out = sys.stdout 
    else:
        if opts.output.endswith(".gz"):
            f_out = gzip.open(opts.output, "wt")
        else:
            f_out = open(opts.output, "w")
        
    # Features
    opts.features = set(filter(bool, map(str.strip, opts.features.split(","))))
    hmmsearch_attribute_label = opts.pyhmmsearch_attribute_label
    if opts.annotation_format in {
            "pykofamsearch", 
            "pykofamsearch-reformatted", 
            "pyhmmsearch",
            "pyhmmsearch-reformatted",
            }:
        hmmsearch_attribute_label = "KOfam"
        
    if opts.annotation_format == "veba":
        id_to_annotation = get_annotations_from_veba(opts.annotations)
    else:
        if opts.annotations.endswith(".gz"):
            f = gzip.open(opts.annotations, "rt")
        else:
            f = open(opts.annotations, "r") 
            
        # Remove header
        if opts.annotation_format in {
            "pykofamsearch", 
            "pykofamsearch-reformatted", 
            "pyhmmsearch",
            "pyhmmsearch-reformatted",
            }:
            f = next(f)
            
        if opts.annotation_format in {
            "pykofamsearch", 
            "pykofamsearch-no-header", 
            "pyhmmsearch",
            "pyhmmsearch-no-header",
            }:
            id_to_annotation = get_annotations_from_pyhmmsearch(f, hmmsearch_attribute_label)
        else:
            id_to_annotation = get_annotations_from_pyhmmsearch_reformatted(f, hmmsearch_attribute_label)


    # Update GFF
    for line in f_in:
        line = line.strip()
        if line.startswith("#"):
            print(line, file=f_out)
        else:
            id_contig, source, feature, start, end, score, strand, frame, attributes_unformatted = line.split("\t")
            if feature in opts.features:
                attributes = list()
                for attribute in attributes_unformatted.split(";"):
                    attribute = attribute.strip()
                    if attribute:
                        k, v = attribute.split("=")
                        attributes.append("{}={}".format(k, v))
                        if k == opts.identifier_attribute:
                            database_to_annotation = id_to_annotation[v]
                            for database, annotation in database_to_annotation.items():
                                if pd.notnull(annotation):
                                    if ";" in annotation:
                                        print("Replacing semicolon with hyphen in annotation: {}".format(annotation), file=sys.stderr)
                                        annotation = annotation.replace(";","-")
                                    if "=" in annotation:
                                        print("Replacing equal sign with hyphen in annotation: {}".format(annotation), file=sys.stderr)
                                        annotation = annotation.replace("=","-")
                                    attributes.append("{}={}".format(database, annotation))
                print("\t".join([id_contig, source, feature, start, end, score, strand, frame, ";".join(attributes)]), file=f_out)
            else:
                print(line, file=f_out)

    # Close files
    if f_in is not sys.stdin:
        f_in.close()
    if f_out is not sys.stdout: 
        f_out.close()


if __name__ == "__main__":
    main()
    
            