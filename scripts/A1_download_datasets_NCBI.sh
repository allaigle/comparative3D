#!/bin/bash
#
####################################################################################################################
#                                    Download genome from NCBI                         
#                        using datasets (NCBI), bash and the SLURM scheduler                                                                             
#
# [USAGE] bash A1_download_datasets_NCBI.sh
#
# [Input] GenBank assembly
# [Output] FASTA
#
# [Adaptation & History]
# Sept,Nov 2024 Alice Laigle (alice.laigle@gmail.com)
# April 2025 AL - made it as array
#
###################################################################################################################
#
#SBATCH --time=168:00:00 # Each task takes max 7 days
#SBATCH --mem-per-cpu=8G # Each task uses max 8GB memory
#SBATCH --cpus-per-task=1 # Each task uses 1 CPU
#SBATCH --array=1-60%10 # submit 60 tasks with task ID 1,2,...,60.

# Variables - TO BE CHANGED
basePATH="/data/alicel/chapter2"
datasets="${basePATH}/4_tools/datasets" # location tool
INFILE="${basePATH}/0_data/0_collection/allDarwin_phylum_species_assembly_SRA.filtered_60.txt"
##### {INFILE} contains tab-separated lines that look like #####
##### basidio	Agentilis	GCA_965113095.1 #####
### this file contains 19 lines with 3 arguments per line
### line <i> contains arguments for run <i>
# get first argument (abbreviated phylum)
PHYLUM=$(cat $INFILE | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')
# get second argument (abbreviated species name)
SPECIES=$(cat $INFILE | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $2}')
# get second argument (GenBank assembly)
ASSEMBLY=$(cat $INFILE | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $3}')

# download genome
mkdir -p ${basePATH}/0_data/1_genomes/slurm
cd ${basePATH}/0_data/1_genomes/slurm

sbatch -J A1_${SPECIES} -c 1 --mem=8MB --wrap="mkdir -p ../${PHYLUM}/${SPECIES} ;\
     $datasets download genome accession $ASSEMBLY --include genome,seq-report \
     --filename ../${PHYLUM}/${SPECIES}/${SPECIES}_${ASSEMBLY}_dataset.zip ;\
     cd ../${PHYLUM}/${SPECIES} ;\
     unzip ${SPECIES}_${ASSEMBLY}_dataset.zip"
