#!/usr/bin/env python

"""
# Original Source: 
# https://github.com/cruizperez/MicrobeAnnotator/blob/master/microbeannotator/pipeline/ko_mapper.py
# Forked on 2023.10.18

If you use this software, please cite the original publication: 
    Ruiz-Perez, C.A., Conrad, R.E. & Konstantinidis, K.T. 
    MicrobeAnnotator: a user-friendly, comprehensive functional annotation pipeline for microbial genomes. 
    BMC Bioinformatics 22, 11 (2021). https://doi.org/10.1186/s12859-020-03940-5

########################################################################
# Author:       Carlos A. Ruiz Perez
# Email:        cruizperez3@gatech.edu
# Intitution:   Georgia Institute of Technology
# Version:      1.0.0
# Date:         Nov 13, 2020

# Description: Maps protein KO information with their respective modules
# and calculates the completeness percentage of each module present.
########################################################################
"""

################################################################################
"""---0.0 Import Modules---"""

import sys, os, re, pickle, argparse, gzip
from collections import OrderedDict, defaultdict
import pandas as pd

__version__ = "2023.10.23"
__program__ = os.path.split(sys.argv[0])[-1]

################################################################################
def load_pickle(path):
    with open(path,"rb") as f:
        obj = pickle.load(f)
    return obj
    
def ko_match(string):
    """ 
    Looks if string has the form K[0-9]{5}
    
    Arguments:
        string {string} -- String to test
    
    Returns:
        [bool] -- String has form or not
    """
    if re.search(r'^K[0-9]{5}$', string) is not None:
        return 1
    else:
        return 0

def split_compound(string, comp_type):
    """[summary]
    
    Arguments:
        string {[type]} -- [description]
        comp_type {[type]} -- [description]
    
    Returns:
        [type] -- [description]
    """
    if comp_type == 'compound':
        return re.split('[-+]', string)
    elif comp_type == 'and_comp':
        return string.split('_')
    elif comp_type == 'or_comp':
        return string.split(',')

def process_compounds(string, comp_type, ko_list):
    string = split_compound(string, comp_type)
    proteins_required = len(string)
    proteins_present = 0
    for option in string:
        if option == '':
            proteins_required -= 1
        elif '+' in option or '-' in option:
            compound = split_compound(option, 'compound')
            proteins_in_compound = len(compound)
            present_in_compound = 0
            for sub_option in compound:
                if ko_match(sub_option) > 0 and sub_option in ko_list:
                    present_in_compound += 1
            proteins_present += present_in_compound/proteins_in_compound
        else:
            if ko_match(option) > 0 and option in ko_list:
                proteins_present += 1
    return proteins_present, proteins_required



def get_kegg_module_information(path):
                                     
    # Get genome and module ids
    module_information = OrderedDict()

    # Get module correspondence
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            id_module, name, group, color = line.split("\t")
            module_information[id_module] = {"module_name":name, "pathway_group":group}

    return pd.DataFrame(module_information).T

def read_ko_list(path):
    ko_list = []

    if path.endswith(".gz"):
        f = gzip.open(path, "rt")
    else:
        f = open(path, "r")
    for line in f:
        ko_list.append(line.strip())
    f.close()
    return set(ko_list)

def regular_module_mapper(ko_list, module_dictionary):#, ko_list_file):
    # ko_list = read_ko_list(ko_list_file)
    regular_module_completenes = []
    for module, final_steps in module_dictionary.items():
        complete_steps = 0
        total_module_steps = len(final_steps)
        for proteins in final_steps.values():
            score_for_step = 0
            steps_per_option = 100
            match = False
            for option in proteins:
                match = False
                # Search for single gene option
                if ko_match(option) > 0:
                    if option in ko_list:
                        score_for_step = 1
                        match = True
                    if match == True:
                        steps_per_option = 1
                    elif 1 < steps_per_option:
                        steps_per_option = 1
                elif '%' in option:
                    option = option.replace(')', '')
                    option = option.split('-%')
                    mandatory = split_compound(option[0], 'compound')
                    proteins_required = len(mandatory) + 1
                    proteins_present = 0
                    for prot in mandatory:
                        if ko_match(prot) > 0 and prot in ko_list:
                            proteins_present += 1
                    for prot in option[1].split(','):
                        if ko_match(prot) > 0 and prot in ko_list:
                            proteins_present += 1
                            break
                    if proteins_present/proteins_required > score_for_step:
                        score_for_step = proteins_present/proteins_required
                        match = True
                    if match == True:
                        steps_per_option = 1
                    elif 1 < steps_per_option:
                        steps_per_option = 1
                elif '_' in option and ',' in option:
                    option = sorted(split_compound(option, 'and_comp'), key=len)
                    highest_score = 0
                    for element in option:
                        if ko_match(element) > 0 and element in ko_list:
                            highest_score += 1
                        elif ',' in element:
                            element = sorted(element.split(","), key=len)
                            for sub_element in element:
                                if ko_match(sub_element) > 0 and sub_element in ko_list:
                                    highest_score += 1
                                    break
                                elif '+' in sub_element:
                                    proteins_present, proteins_required = process_compounds(sub_element, 
                                    'compound', ko_list)
                                    highest_score += proteins_present/proteins_required
                    if highest_score > score_for_step:
                        score_for_step = highest_score
                        match = True
                    if match == True:
                        steps_per_option = len(option)
                    elif len(option) < steps_per_option:
                        steps_per_option = len(option)
                elif len(option.split(",")) > 1:
                    match = False
                    option = sorted(option.split(","), key=len)
                    highest_score = 0
                    steps_to_add = 0
                    for sub_option in option:
                        if ko_match(sub_option) > 0 and sub_option in ko_list:
                            if 1 > highest_score:    
                                highest_score = 1
                                steps_to_add = 1
                        elif '+' in sub_option or '-' in sub_option:
                            proteins_present, proteins_required = process_compounds(sub_option, 'compound', ko_list)
                            if proteins_present/proteins_required > highest_score:
                                highest_score = proteins_present/proteins_required
                                steps_to_add = 1
                    if highest_score > score_for_step:
                        score_for_step = highest_score
                        match = True
                    if match == True:
                        steps_per_option = steps_to_add
                    elif steps_to_add< steps_per_option:
                        steps_per_option = steps_to_add
                elif '_' in option:
                    match = False
                    proteins_present, proteins_required = process_compounds(option, 'and_comp', ko_list)
                    if proteins_present/proteins_required > score_for_step:
                        score_for_step = proteins_present
                        match = True
                    if match == True:
                        steps_per_option = proteins_required
                    elif proteins_required < steps_per_option:
                        steps_per_option = proteins_required
                elif "+" in option or "-" in option:
                    match = False
                    highest_score = 0
                    proteins_present, proteins_required = process_compounds(option, 'compound', ko_list)
                    if proteins_present/proteins_required > score_for_step:
                        score_for_step = proteins_present/proteins_required
                        match = True
                    if match == True:
                        steps_per_option = 1
                    elif proteins_required < steps_per_option:
                        steps_per_option = 1
                else:
                    print("Unrecognized module {}. Check your database.".format(option), file=sys.stderr)
            complete_steps += score_for_step
            if steps_per_option > 50:
                steps_per_option = 1
            if steps_per_option > 1:
                total_module_steps += steps_per_option - 1
        regular_module_completenes.append((module, (complete_steps/total_module_steps)))
    return regular_module_completenes

def bifurcating_module_mapper(ko_list, module_dictionary):#, ko_list_file):
    # ko_list = read_ko_list(ko_list_file)
    bifurcating_module_completenes = []
    for module, versions in module_dictionary.items():
        module_highest = 0
        for version, total_steps in versions.items():
            completed_steps = 0
            total_version_steps = len(total_steps)
            for proteins in total_steps.values():
                score_for_step = 0
                steps_per_option = 100
                match = False
                if isinstance(proteins, (list)):
                    protein = sorted(proteins, key=len)
                    for option in protein:
                        match = False
                        if ko_match(option) > 0:
                            if option in ko_list:
                                score_for_step = 1
                                match = True
                            if match == True:
                                steps_per_option = 1
                            elif 1 < steps_per_option:
                                steps_per_option = 1
                        elif '_' in option and ',' in option:
                            option = sorted(split_compound(option, 'and_comp'), key=len)
                            highest_score = 0
                            for element in option:
                                if ko_match(element) > 0 and element in ko_list:
                                    highest_score += 1
                                elif ',' in element:
                                    element = sorted(element.split(","), key=len)
                                    score_sub_option = 0
                                    for sub_element in element:
                                        if ko_match(sub_element) > 0 and sub_element in ko_list:
                                            highest_score += 1
                                            break
                                        elif '+' in sub_element:
                                            proteins_present, proteins_required = process_compounds(sub_element, 
                                            'compound', ko_list)
                                            highest_score += proteins_present/proteins_required
                            if highest_score > score_for_step:
                                score_for_step = highest_score
                                match = True
                            if match == True:
                                steps_per_option = len(option)
                        elif '_' in option:
                            match = False
                            proteins_present, proteins_required = process_compounds(option, 'and_comp', ko_list)
                            if proteins_present/proteins_required > score_for_step:
                                score_for_step = proteins_present
                                match = True
                            if match == True:
                                steps_per_option = proteins_required
                            elif proteins_required < steps_per_option:
                                steps_per_option = proteins_required
                        elif '+' in option or '-' in option:
                            match = False
                            highest_score = 0
                            proteins_present, proteins_required = process_compounds(option, 'compound', ko_list)
                            if proteins_present/proteins_required > score_for_step:
                                score_for_step = proteins_present/proteins_required
                                match = True
                            if match == True:
                                steps_per_option = 1
                            elif proteins_required < steps_per_option:
                                steps_per_option = 1
                elif ko_match(proteins) > 0:
                    match = False
                    if proteins in ko_list:
                        score_for_step = 1
                        match = True
                    if match == True:
                        steps_per_option = 1
                    elif 1 < steps_per_option:
                        steps_per_option = 1
                elif ',' in proteins:
                    options = split_compound(proteins, 'or_comp')
                    for option in options:
                        if ko_match(option) > 0 and option in ko_list:
                            score_for_step = 1
                            match = True
                        if match == True:
                            steps_per_option = 1
                        elif 1 < steps_per_option:
                            steps_per_option = 1
                elif '+' in proteins or '-' in proteins:
                    match = False
                    highest_score = 0
                    proteins_present, proteins_required = process_compounds(proteins, 'compound', ko_list)
                    if proteins_present/proteins_required > score_for_step:
                        score_for_step = proteins_present/proteins_required
                        match = True
                    if match == True:
                        steps_per_option = 1
                    elif proteins_required < steps_per_option:
                        steps_per_option = 1
                else:
                    print("Unreccognized module {}. Check your database.".format(proteins), file=sys.stderr)
                completed_steps += score_for_step
                if steps_per_option > 50:
                    steps_per_option = 1
                if steps_per_option > 1:
                    total_version_steps += steps_per_option - 1
            if completed_steps/total_version_steps > module_highest:
                module_highest = completed_steps/total_version_steps
        bifurcating_module_completenes.append((module, module_highest))
    return bifurcating_module_completenes

def structural_module_mapper(ko_list, module_dictionary):#, ko_list_file):
    # ko_list = read_ko_list(ko_list_file)
    structural_module_completeness = []
    for module, components in module_dictionary.items():
        score_for_components = 0
        module_proteins_present = 0
        module_proteins_required = 0
        for proteins in components:
            if isinstance(proteins, (list)):
                highest_score = 0
                proteins_to_add = 0
                steps_to_add = 100
                for option in proteins:
                    if '_' in option and ',' in option:
                        option = sorted(split_compound(option, 'and_comp'), key=len)
                        proteins_present_option = 0
                        proteins_required_option = 0
                        for element in option:
                            if ko_match(element) > 0:
                                if element in ko_list:
                                    proteins_present_option += 1
                                    proteins_required_option += 1
                                else:
                                    proteins_required_option += 1
                            elif ',' in element:
                                element = sorted(element.split(","), key=len)
                                score_sub_element = 0
                                proteins_present_sub_element = 0
                                proteins_required_sub_element = 0
                                for sub_element in element:
                                    if ko_match(sub_element) > 0:
                                        if sub_element in ko_list:
                                            score_sub_element = 1
                                            proteins_present_sub_element += 1
                                            proteins_required_sub_element += 1
                                        elif 1 < proteins_required_sub_element and score_sub_element == 0:
                                            proteins_required_sub_element = 1
                                    elif '+' in sub_element or '-' in sub_element:
                                        proteins_present, proteins_required = process_compounds(sub_element, 
                                        'compound', ko_list)
                                        if proteins_present/proteins_required > score_sub_element:
                                            score_sub_element = proteins_present/proteins_required
                                            proteins_present_sub_element = proteins_present
                                            proteins_required_sub_element = proteins_required
                                        elif proteins_required < proteins_required_sub_element and score_sub_element == 0:
                                            proteins_required_sub_element = proteins_required
                                proteins_present_option += proteins_present_sub_element
                                proteins_required_option += proteins_required_sub_element
                        if proteins_present_option/proteins_required_option > highest_score:
                            highest_score = proteins_present_option/proteins_required_option
                            proteins_to_add = proteins_present_option
                            steps_to_add = proteins_required_option
                        elif proteins_required_option < steps_to_add and highest_score == 0:
                            steps_to_add = proteins_required_option
                    elif '+' in option or '-' in option:
                        proteins_present, proteins_required = process_compounds(option, 'compound', ko_list)
                        if proteins_present/proteins_required > highest_score:
                            highest_score = proteins_present/proteins_required
                            steps_to_add = proteins_required
                            proteins_to_add = proteins_present
                        elif proteins_required < steps_to_add and highest_score == 0:
                            steps_to_add = proteins_required
                module_proteins_present += proteins_to_add
                module_proteins_required += steps_to_add
            elif ko_match(proteins) > 0:
                if proteins in ko_list:
                        module_proteins_present += 1
                        module_proteins_required += 1
                else:
                    module_proteins_required += 1
            elif ',' in proteins:
                proteins = sorted(proteins.split(","), key=len)
                score = 0
                proteins_present_option = 0
                proteins_required_option = 100
                for element in proteins:
                    if ko_match(element) > 0:
                        if element in ko_list:
                            proteins_required_option = 1
                            proteins_present_option = 1
                            break
                        elif 1 < proteins_required_option and score == 0:
                            proteins_required_option = 1
                    elif '+' in element or '-' in element:
                        proteins_present, proteins_required = process_compounds(element, 'compound', ko_list)
                        if proteins_present/proteins_required > score:
                            score = proteins_present/proteins_required
                            proteins_present_option = proteins_present
                            proteins_required_option = proteins_required
                        elif proteins_required < proteins_required_option and score == 0:
                            proteins_required_option = proteins_required
                module_proteins_present += proteins_present_option
                module_proteins_required += proteins_required_option
            elif '+' in proteins or '-' in proteins:
                proteins_present, proteins_required = process_compounds(proteins, 'compound', ko_list)
                module_proteins_present += proteins_present
                module_proteins_required += proteins_required  
            else:
                print("Unreccognized module {}. Check your database.".format(proteins), file=sys.stderr)
        score_for_components = module_proteins_present/module_proteins_required
        structural_module_completeness.append((module, score_for_components))
    return structural_module_completeness



################################################################################
"""---2.0 Main Function---"""

def main():
    # Setup parser for arguments.

    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i ko_table.tsv [-k <GENOME_1.ko_list [GENOME_2.ko_list ... GENOME_N.ko_list]] > -d <database_directory -o <output.tsv[.gz]> -x ko_list".format(__program__)

    epilog = "Josh L. Espinoza's fork from MicrobeAnnotator. Please cite the following: https://doi.org/10.1186/s12859-020-03940-5"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('-i', '--ko_table', help='path/to/ko_table.tsv following [id_genome]<tab>[id_ko], No header.  Cannot be used with --ko_lists')
    parser.add_argument('-k', '--ko_lists', nargs='+', help='Space-delimited list of filepaths where each file represents a genome and each line in the file is a KO id.  Cannot be used with --ko_table')
    parser.add_argument('-o', '--output',  default="stdout", help='Output file for module completion ratios [Default: stdout]')
    parser.add_argument("-d", '--database_directory', required=True, help='path/to/database_directory with pickle files')
    parser.add_argument("-x", '--extension', default = "ko_list", help='File extension for --ko_lists [Default: ko_list (e.g., MAG_1.ko_list)]')
    parser.add_argument("-t", '--transpose', action='store_true', help='Transpose output format.  Default rows=MCR, columns=Genomes.')
    parser.add_argument("-m", '--no_missing_modules', action='store_true', help='Do not include modules that are missing from all genomes')
    parser.add_argument('--no_module_names', action='store_true', help='Do not include module names in output')
    parser.add_argument('--no_pathway_groups', action='store_true', help='Do not include pathway group names in output')


    opts = parser.parse_args()

    assert bool(opts.ko_table) != bool(opts.ko_lists), "Must provide KOs as either a tsv table (--ko_table) or a list of KO ids in different files (--ko_lists)"
    if opts.ko_table == "stdin":
        opts.ko_table = sys.stdin 

    if opts.output == "stdout":
        opts.output = sys.stdout

    # ----------------------------
    # Import all modules from dictionaries
    regular_modules = load_pickle(os.path.join(opts.database_directory, "KEGG_Regular_Module_Information.pkl"))
    bifurcating_modules = load_pickle(os.path.join(opts.database_directory, "KEGG_Bifurcating_Module_Information.pkl"))
    structural_modules = load_pickle(os.path.join(opts.database_directory, "KEGG_Structural_Module_Information.pkl"))
             
    # Get KEGG ortholog lists
    if opts.ko_table:
        genome_to_kos = defaultdict(set)
        df = pd.read_csv(opts.ko_table, sep="\t", index_col=None, header=None).iloc[:,:2]
        for _, (id_genome, id_ko) in df.iterrows():
            id_genome = str(id_genome)
            genome_to_kos[id_genome].add(id_ko)
    else:
        genome_to_kos = OrderedDict()
        for path in opts.ko_lists:
            id_genome = path[:-(len(opts.extension) + 1)]
            id_genome = str(id_genome)
            genome_to_kos[id_genome] = read_ko_list(path)

    # Get module information
    df_module_information = get_kegg_module_information(
        path=os.path.join(opts.database_directory, "KEGG_Module_Information.txt"),
        )

    # Calculate module completion ratios
    module_completion_ratios = OrderedDict()
    for id_genome, ko_list in genome_to_kos.items():
        regular_completeness = regular_module_mapper(ko_list, regular_modules)
        bifurcating_completeness = bifurcating_module_mapper(ko_list, bifurcating_modules)
        structural_completeness = structural_module_mapper(ko_list, structural_modules)
        final_completeness = regular_completeness + bifurcating_completeness + structural_completeness
        module_completion_ratios[id_genome] = pd.Series(dict(final_completeness))
    df_mcr = pd.DataFrame(module_completion_ratios)

    # Format output
    df_output = pd.concat([df_module_information, df_mcr], axis=1)
    if opts.no_missing_modules:
        df_output = df_output.loc[df_mcr.sum(axis=1)[lambda x: x > 0].index]
    df_output.index.name = "id_kegg-ortholog"

    if opts.no_module_names:
        df_output = df_output.drop(["module_name"], axis=1)
    if opts.no_pathway_groups:
        df_output = df_output.drop(["pathway_group"], axis=1)
    else:
        df_output = df_output.reset_index(drop=False).set_index(["pathway_group", "id_kegg-ortholog"])
        
    df_output.to_csv(opts.output, sep="\t")

if __name__ == "__main__":
    main()
