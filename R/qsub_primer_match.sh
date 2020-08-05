#!/bin/bash
# Commands starting with  '#$' are interpreted by SGE
# Shell to be used for the job
#$ -S /bin/bash
# User to be informed
#$ -M vaulot@gmail.com
# Export all environment variable
#$ -V
# Send a message by email  at beginning (b), end (e) and abort (a) of job
#$ -m bea
# Standard output.  Can use '-j y' to add stderr with stdout
#$ -o repl
# Send the commande from the curent directory where the script reside
#$ -cwd
## Queue settings (can be set at runtime in the qsub submition line)
#$ -q short.q
## [Optional] to ask between 1 and 4 processors
##$ -pe thread 1-16
## [Recommended] Ask for at least 1GB of memory for the job
##$ -l mem_free=1G
## [Recommended] Kill the job if it eats 4GB per slot (i.e. thread)
## here it means between 4GB*1 thread=4GB and 4GB*4 threads=16GB
## $ -l h_vmem=4G - Desactivate


# Submitted with parameter providing the primer_set_id
# qsub qsub_primer_match.sh 36


# Print information about the current job
# ---------------------------------------
# Print beginning of job
echo "$(date) job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME" 1>&2
# Print the number of threads used
echo "NSLOTS = $NSLOTS" 1>&2


cd /projet/umr7144/dipo/vaulot/pr2/primers
/opt/6.x/R-3.5.1/bin/Rscript --no-save --no-restore script_primers_pr2_match.R -s $1 > script_primers_pr2_match_set_$1.out

# Get summary of job stats
# ------------------------
# Get (max) memory consumption
echo "$(qstat -f -j $JOB_ID | grep vmem)" 1>&2
# Print end of job. If not printed, job died before the end.
echo "$(date) job $JOB_NAME done" 1>&2