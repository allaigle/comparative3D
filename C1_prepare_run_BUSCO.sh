#!/bin/bash
#SBATCH -c 1 # processor' number
#SBATCH --mem 2MB # memory

# Adaptation Alice Laigle from Toby Baril's script.
# [USAGE] sbatch -J C1 -c 1 --mem=2MB C1_prepare_run_BUSCO.sh 

# Copy entire FASTA genomes (darwin subset + 3 zoopagomycota; not filtered wo500)

basePATH="/data/alicel/chapter2"
workdir="${basePATH}/1_dataTreatment/2_phyloTree"
INFILE="${basePATH}/0_data/0_collection/phylotree_phylum_species_assembly.63.txt"

cd $workdir 

while read PHYLUM SPECIES ASSEMBLY _ ; do

  cp ${basePATH}/0_data/1_genomes/${PHYLUM}/${SPECIES}/ncbi_dataset/data/${ASSEMBLY}/${ASSEMBLY}*.fna \
    ${workdir}/0_inputs/genomes/${SPECIES}.${ASSEMBLY}.fna  

done < "$INFILE"

mkdir -p ${workdir}/1_busco && cd ${workdir}/1_busco

ln -s /data/alicel/chapter2/3_scripts/C2_run_busco_TBaril.sh C2_run_busco_TBaril.sh

# Created buscoParameters.63.txt 
##### buscoParameters.txt contains tab-separated lines that look like #####
##### Abisporus.GCA_943193715.1.fna   BUSCO_Abisporus #####
### this file contains 63 lines with 2 arguments per line

conda activate BUSCO

bash C2_run_busco_TBaril.sh

mkdir -p ${workdir}/2_busco/slurm ; mv slurm-* slurm/