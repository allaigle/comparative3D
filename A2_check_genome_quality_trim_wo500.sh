#!/bin/bash
#
####################################################################################################################
#                          Check genome quality and trim under 500kb                  
#                   using samtools, seqtk, bash and the SLURM scheduler                                                                             
#
# [GAOL] Create genome indices and remove all chromosomes, SUPER contigs or scaffolds under 500kb.
#
# [USAGE] sbatch A2_check_genome_quality_trim_wo500.sh
#
# [Adaptation & History]
# Nov 2024 Alice Laigle (alice.laigle@unine.ch)
# Apr 2025 AL - Remove arguments and add INFILE for while loop
# 
# NOTES:
#    1. Genomes are already downloaded from NCBI. If needed, run 'A1_download_datasets_NCBI.sh'
#    2. samtools, seqtk are already installed.  
###################################################################################################################
#
#SBATCH -c 1 # processor' number
#SBATCH --mem 2MB # memory

# Variables - TO BE CHANGED
basePATH="/data/alicel/chapter2"
workdir="${basePATH}/0_data/1_genomes/wo500"
INFILE="${basePATH}/0_data/0_collection/allDarwin_phylum_species_assembly_SRA.filtered_60.txt"
export PATH=$PATH:/data/alicel/chapter2/4_tools/seqtk

mkdir -p ${workdir} && cd ${workdir} 

while read PHYLUM SPECIES ASSEMBLY _ ; do

  # Input - FASTA
  ln -sf ${basePATH}/0_data/1_genomes/${PHYLUM}/${SPECIES}/ncbi_dataset/data/${ASSEMBLY}/${ASSEMBLY}*.f*a \
    ${workdir}/${SPECIES}_${ASSEMBLY}.fasta

  # Create genome index
  samtools faidx ${workdir}/${SPECIES}_${ASSEMBLY}.fasta 

  # Exclude every sequence under 500kb
  seqtk seq -L 500000 ${workdir}/${SPECIES}_${ASSEMBLY}.fasta \
    > ${workdir}/${SPECIES}_${ASSEMBLY}.wo500.fasta

  # Create the 2nd index
  samtools faidx ${workdir}/${SPECIES}_${ASSEMBLY}.wo500.fasta

  # Check the number of chr/scaffold/contigs between both indices
  wc -l ${workdir}/${SPECIES}_${ASSEMBLY}.fasta.fai \
    ${workdir}/${SPECIES}_${ASSEMBLY}.wo500.fasta.fai 

done < $INFILE

# Give accesses
chmod -R 775 ${workdir}