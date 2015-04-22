PATH=$PATH:$HOME/bin
module load slurm
module load python/anaconda-2.1.0 

sbatch ./scripts/_runScript.sh "$@"
