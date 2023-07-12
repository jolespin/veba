# __version__ = "2023.7.11"

VEBA_VERSION=$(head -n 1 ../VERSION)
VDB_VERSION=$(head -n 2 ../VERSION | tail -n 1)

>SCRIPT_VERSIONS
echo "VEBA __version__ = \"${VEBA_VERSION}\"" >> SCRIPT_VERSIONS
echo "VEBA_DATABASE __version__ = \"${VDB_VERSION}\"" >> SCRIPT_VERSIONS

for FP in *.py scripts/*.py scripts/*.r ../install/*.sh ../install/docker/*.sh; do V=$(grep "__version__ =" ${FP}); echo "${FP} ${V}"; done >> SCRIPT_VERSIONS
