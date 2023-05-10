N=install
N_JOBS=1
mkdir -p logs
CMD="bash install_veba.sh"
sbatch -J ${N} -p ind-shared -N 1 -c ${N_JOBS} --ntasks-per-node=1 -A jcl119 -o logs/${N}.o -e logs/${N}.e --export=ALL -t 12:00:00 --mem=15G --wrap="${CMD}"
