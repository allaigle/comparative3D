#!/bin/bash

#SBATCH --time=168:00:00 # Each task takes max 7 days
#SBATCH --mem-per-cpu=8G # Each task uses max 8GB memory
#SBATCH --cpus-per-task=2 # Each task uses 2 CPUs
#SBATCH --array=1-758%10 # submit 10 tasks with task ID 1,2,...,758.

# Adaptation Alice Laigle from Toby Baril's script.
# [USAGE] bash C6_run_mafft_TBaril.sh

# set directories
data_dir=/data/alicel/chapter2/1_dataTreatment/2_phyloTree/2_busco_aa

# get parameters - this is a file with the name of each AA fasta file
param_store=/data/alicel/chapter2/1_dataTreatment/2_phyloTree/alignmentParameters.txt
##### param_store contains single column containing the name for each file #####
##### gene1_aa.fasta #####
### this file contains X lines with 1 arguments per line
### line <i> contains arguments for run <i>
# get first argument (name of alignment file)
arg_one=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')

# run mafft
srun mafft --auto ${data_dir}/${arg_one} > /data/alicel/chapter2/1_dataTreatment/2_phyloTree/3_alignments/${arg_one%.fasta}.aln
