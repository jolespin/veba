# MicrobeAnnotator-KEGG

**If this is used in any way, please cite the source publication:** 

Ruiz-Perez, C.A., Conrad, R.E. & Konstantinidis, K.T. MicrobeAnnotator: a user-friendly, comprehensive functional annotation pipeline for microbial genomes. BMC Bioinformatics 22, 11 (2021). https://doi.org/10.1186/s12859-020-03940-5

**This data has been incorporated from the following source:** 

https://github.com/cruizperez/MicrobeAnnotator/tree/master/microbeannotator/data

**File Descriptions:**

* `KEGG_Regular_Module_Information.pkl` - Python dictionary of regular modules from `MicrobeAnnotator` of `{id_module:structured_kegg_orthologs}`
* `KEGG_Bifurcating_Module_Information.pkl` - Python dictionary of bifurcating modules from `MicrobeAnnotator` of `{id_module:structured_kegg_orthologs}`
* `KEGG_Structural_Module_Information.pkl` - Python dictionary of structural modules from `MicrobeAnnotator` of `{id_module:structured_kegg_orthologs}`
* `KEGG_Module_Information.txt` - - Table containing KEGG ortholog, higher level categories, and module color
* `KEGG_Module-KOs.pkl` - Flattened dictionary which includes `{id_module:{KO_1, KO_2, ..., KO_M}`.  Note: This is not structured and should be used cautiously as KEGG modules and completion calculations are complex.  Generated with the Python code below:

```python
import pickle, glob, os

kegg_directory = "{}/MicrobeAnnotator_KEGG/".format(os.environ["VEBA_DATABASE"]

delimiters = [",","_","-","+"]

# Load MicrobeAnnotator KEGG dictionaries
module_to_kos__unprocessed = defaultdict(set)
for fp in glob.glob(os.path.join(kegg_directory, "*.pkl")):
    with open(fp, "rb") as f:
        d = pickle.load(f)
        
    for id_module, v1 in d.items():
        if isinstance(v1, list):
            try:
                module_to_kos__unprocessed[id_module].update(v1)
            except TypeError:
                for v2 in v1:
                    module_to_kos__unprocessed[id_module].update(v2)
        else:
            for k2, v2 in v1.items():
                if isinstance(v2, list):
                    try:
                        module_to_kos__unprocessed[id_module].update(v2)
                    except TypeError:
                        for v3 in v2:
                            module_to_kos__unprocessed[id_module].update(v3)

# Flatten the KEGG orthologs
module_to_kos__processed = dict()
for id_module, kos_unprocessed in module_to_kos__unprocessed.items():
    kos_processed = set()
    for id_ko in kos:
        composite=False
        for sep in delimiters:
            if sep in id_ko:
                id_ko = id_ko.replace(sep,";")
                composite = True
        if composite:
            kos_composite = set(map(str.strip, filter(bool, id_ko.split(";"))))
            kos_processed.update(kos_composite)
        else:
            kos_processed.add(id_ko)
    module_to_kos__processed[id_module] = kos_processed


# Write
with open(os.path.join(kegg_directory, "KEGG_Module-KOs.pkl"), "wb") as f:
    pickle.dump(module_to_kos__processed, f)
```