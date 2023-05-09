N="database_config"
N_JOBS=1
DATABASE_DIRECTORY=/expanse/projects/jcl110/db/veba/VDB_v4-test
CMD="source activate test-VEBA-database_env && bash download_databases-lite.sh ${DATABASE_DIRECTORY}"
sbatch -J ${N} -p ind-shared -N 1 -c ${N_JOBS} --ntasks-per-node=1 -A jcl110 -o logs/${N}.o -e logs/${N}.e --export=ALL -t 14:00:00 --mem=16G --wrap="${CMD}"

# CMD="source activate test-VEBA-database_env && bash download_databases.sh ${DATABASE_DIRECTORY}"
# sbatch -J ${N} -p ind-shared -N 1 -c ${N_JOBS} --ntasks-per-node=1 -A jcl110 -o logs/${N}.o -e logs/${N}.e --export=ALL -t 14:00:00 --mem=64G --wrap="${CMD}"
