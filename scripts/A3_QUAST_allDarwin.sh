#!/bin/bash
#
####################################################################################################################
#                              Check genome quality                  
#                   using QUAST, bash and the SLURM scheduler                                                                             
#
# [USAGE] sbatch -J A3 -c 1 --mem=16MB A3_QUAST_allDarwin.sh
#
# [Input] FASTA
# [Output] reports
#
# [Adaptation & History]
# Nov 2024 Alice Laigle (alice.laigle@gmail.com)
# Apr 2025 AL - made for all at once, not through species loop
# 
# NOTE: QUAST is already installed in its homonym conda env.  
###################################################################################################################
#
#SBATCH -c 1 # processor' number
#SBATCH --mem 16MB # memory

# Variables - TO BE CHANGED
basePATH="/data/alicel/chapter2"
workdir="${basePATH}/1_dataTreatment/0_QUAST"
QUAST="${basePATH}/4_tools/quast/quast.py"
INFILE="${basePATH}/0_data/0_collection/allDarwin_phylum_species_assembly_SRA.filtered_60.txt"

cd ${workdir}
mkdir -p ${workdir}/fasta

# INPUTS
while read _ SPECIES ASSEMBLY _; do

  ln -sf ${basePATH}/0_data/1_genomes/wo500/${SPECIES}_${ASSEMBLY}.fasta \
    ${workdir}/fasta/${SPECIES}_${ASSEMBLY}.fasta 

  ln -sf ${basePATH}/0_data/1_genomes/wo500/${SPECIES}_${ASSEMBLY}.wo500.fasta\
    ${workdir}/fasta/${SPECIES}_${ASSEMBLY}.wo500.fasta

done < $INFILE

# Run QUAST for all the genomes

conda activate QUAST

$QUAST ${workdir}/fasta/*.fasta --fungus -o allDarwin

conda deactivate

chmod -R 775 ${workdir}
