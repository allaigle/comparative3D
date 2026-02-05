#!/bin/bash

# Adaptation Alice Laigle from Toby Baril's script.
# [USAGE] bash C2_run_busco_TBaril.sh

#SBATCH --time=168:00:00 # Each task takes max 7 days
#SBATCH --mem-per-cpu=8G # Each task uses max 8GB memory
#SBATCH --cpus-per-task=2 # Each task uses 2 CPUs
#SBATCH --array=1-63%10 # submit 63 tasks with task ID 1,2,...,63.

# set directories
data_dir=/data/alicel/chapter2/1_dataTreatment/2_phyloTree/0_inputs/genomes

# get parameters - this is a file with the name of each genome file, and the 
## name of the resultant output file (WITHOUT PATHS)
param_store=/data/alicel/chapter2/1_dataTreatment/2_phyloTree/1_busco/buscoParameters.63.txt
##### param_store contains tab-separated lines that look like #####
##### agabis3.fna	BUSCO_agabis3 #####
### this file contains 63 lines with 2 arguments per line
### line <i> contains arguments for run <i>
# get first argument (name of genome file)
arg_one=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')
# get second argument (name for output directory)
arg_two=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $2}')

# run busco
srun busco -i ${data_dir}/${arg_one} -m genome -l fungi_odb10 -c 2 -o ${arg_two}
