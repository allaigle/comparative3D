#!/bin/bash
#
####################################################################################################################
#                                        Run 3DGB                         
#                            using bash and the SLURM scheduler                                                                             
#
# [USAGE] bash B2_run_3DGB.sh PHYLUM SPECIES SRA
#
# [Input] .fasta, fastq.gz and YML config file
# [Output] check 3DGB github page
#
# [Adaptation & History]
# Sept-Oct 2024 Alice Laigle (alice.laigle@gmail.com)
# 
# NOTE: Includes 'B3_clean_3DGB_workdir.sh' script in the second SBATCH command.
###################################################################################################################
#
# Arguments
PHYLUM=$1 # e.g., PHYLUM="muco"
SPECIES=$2 # e.g., SPECIES="Eparvispora"
SRA=$3 # e.g., SRA="ERR11577535"

# Variable
basePATH="/data/alicel/chapter2"
workdir="${basePATH}/1_dataTreatment/1_3DGB"

#### Environment
conda activate 3DGB
export PATH=/data/alicel/miniconda3/bin:$PATH
# to get "iced" pckg, issues arose when not having these 2 more steps
conda deactivate
conda activate 3DGB 

#### Run 3DGB
cd ${workdir}

sbatch -J 3DGB_${SPECIES} -c 4 --mem=16GB \
  --wrap="snakemake --profile smk_profile_debug -j 4 \
    --configfile ${workdir}/configYML/${PHYLUM}/${SPECIES}_config.yml"

# Cleam 3DGB SPECIES folder
cd ${basePATH}/3_scripts 
sbatch -J clean_3DGB_${SPECIES} -c 1 --mem=2MB \
  --dependency=$(squeue --noheader --format %i --name 3DGB_${SPECIES}) \
  --wrap="bash ${basePATH}/3_scripts/B3_clean_3DGB_workdir.sh \
    ${PHYLUM} ${SPECIES} ${SRA}"
