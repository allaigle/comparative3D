#!/bin/bash
#
####################################################################################################################
#                                Annotatation of TEs                
#                   using EarlGrey, bash and the SLURM scheduler                                                                             
#
# [GAOL] Annotate TEs for each species using the single library created and manually curated.
#
# [USAGE] sbatch -J D6 D6_annotate_TEs_earlGreyAnnotationOnly.sh
#
# [Input] .wo500.fasta & single library
# [Output] TE library, Summary Figures, and TE Quantifications
#
# [Adaptation & History]
# Jan 2025 Alice Laigle (alice.laigle@unine.ch)
# Jun 2025 AL - adapted for all Darwin species
# 
###################################################################################################################
#
#SBATCH --time=168:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=4 
#SBATCH --array=1-60%20


source ~/.bashrc # MIGHT BE CHANGED

# Variables - TO BE CHANGED
basePATH="/data/alicel/chapter2"
workdir="${basePATH}/1_dataTreatment/4_annotations_EarlGrey"
single_library=${workdir}/2_single_library/reduced_library.curated.fa
#export PATH=$PATH:/data/alicel/miniconda3/pkgs/earlgrey-5.0.0-h4ac6f70_1/bin

cd ${workdir} 
mkdir -p ${workdir}/3_final_annotation

# get parameters
param_store=${basePATH}/0_data/0_collection/allDarwin_phylum_species_assembly_SRA.filtered_60.txt
##### param_store contains tab-separated lines that look like #####
##### Abisporus GCA_943193715.1 #####
### this file contains 39 lines with 2 arguments per line
# get first argument (short name of species)
arg_species=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $2}')
# get second argument (name of assembly)
arg_assembly=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $3}')


# fasta input
ln -sf ${basePATH}/0_data/1_genomes/wo500/${arg_species}_${arg_assembly}.wo500.fasta \
  ${workdir}/3_final_annotation/${arg_species}.wo500.fasta
genome="${workdir}/3_final_annotation/${arg_species}.wo500.fasta"

conda activate earlgrey 

earlGreyAnnotationOnly -g $genome -s $arg_species \
  -o ${workdir}/3_final_annotation \
  -l $single_library -t 4 -m yes -e yes