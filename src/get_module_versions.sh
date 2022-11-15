for FP in *.py scripts/*.py; do V=$(grep "__version__ =" ${FP}); echo "${FP} ${V}"; done > ../module_versions
