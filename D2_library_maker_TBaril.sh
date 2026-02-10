#!/bin/bash
#SBATCH --time=168:00:00 # Each task takes max 7 days
#SBATCH --cpus-per-task=16 # Each task uses 16 CPUs
#SBATCH --array=1-8%8 # submit 8 tasks with task ID 1,2,...,8.

# Made and run by Toby Baril
# His EarlGrey github:https://github.com/TobyBaril/EarlGrey

# set directories
data_dir=/data/toby/genomes_alice/1_genomes

# get parameters - this is a file with the name of each genome file, and the 
## name of the resultant output file (WITHOUT PATHS)
param_store=/data/toby/genomes_alice/1_genomes/param_store.txt
### this file contains 8 lines with 1 arguments per line
### line <i> contains arguments for run <i>
# get first argument
arg_one=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')
# get second argument
arg_two=$(echo $arg_one | sed 's/.1_.*//g')

# run earl grey subset
srun earlGreyLibConstruct -g ${data_dir}/${arg_one} -s ${arg_two} -t 16 -o /data/toby/genomes_alice/2_earlGrey
