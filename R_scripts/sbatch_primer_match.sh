#!/bin/bash
#SBATCH -p fast                      # partition
#SBATCH -N 1                         # nombre de nœuds
#SBATCH -n 1                         # nombre de cœurs
#SBATCH --cpus-per-task 1
#SBATCH --mem 32GB                    # mémoire vive pour l'ensemble des cœurs
##SBATCH -t 8-0:00                     # durée maximum du travail (D-HH:MM)
#SBATCH -o slurm.%N.%j.out           # STDOUT
#SBATCH -e slurm.%N.%j.err           # STDERR
#SBATCH --mail-user=vaulot@gmail.com
#SBATCH --mail-type=BEGIN,END


# Partition can be also fast, long, bigmem

## --ntasks=# : Number of "tasks" (use with distributed parallelism).
## --ntasks-per-node=# : Number of "tasks" per node (use with distributed parallelism).
## --cpus-per-task=# : Number of CPUs allocated to each task (use with shared memory parallelism).


# Submitted with parameter providing the primer_set_id
#  sbatch sbatch_primer_match.sh 36



DATE=`date +%Y-%m-%d_%H-%M`

module load r/3.5.1

cd /projet/umr7144/dipo/vaulot/pr2/primers

Rscript --no-save --no-restore script_primers_pr2_match.R -s $1 > script_primers_pr2_match_set_$1.out

