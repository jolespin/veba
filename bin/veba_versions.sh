# __version__ = "2024.1.23"

VEBA_VERSION=$(head -n 1 ../VERSION)
VDB_VERSION=$(head -n 2 ../VERSION | tail -n 1)

>VEBA_SCRIPT_VERSIONS
echo "VEBA __version__ = \"${VEBA_VERSION}\"" >> VEBA_SCRIPT_VERSIONS
echo "VEBA_DATABASE __version__ = \"${VDB_VERSION}\"" >> VEBA_SCRIPT_VERSIONS

for FP in *.py scripts/*.py scripts/*.r ../install/*.sh ../install/docker/*.sh; do V=$(grep "__version__ =" ${FP}); echo "${FP} ${V}"; done >> VEBA_SCRIPT_VERSIONS
