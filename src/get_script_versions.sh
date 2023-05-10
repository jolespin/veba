VEBA_VERSION=$(cat ../VERSION)
echo "VEBA __version__ ${VEBA_VERSION}" > SCRIPT_VERSIONS
for FP in *.py scripts/*.py scripts/*.r; do V=$(grep "__version__ =" ${FP}); echo "${FP} ${V}"; done >> SCRIPT_VERSIONS
