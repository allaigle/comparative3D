#!/bin/bash

#SBATCH -c 1 # processor' number
#SBATCH --mem 2MB # memory

# Adaptation Alice Laigle from Toby Baril's script.
# [USAGE] sbatch -J D7 -c 1 --mem=2MB D7_extract_sequences_to_check.sh

# The goal is to extract sequences for those I doubt about the content
## Cf the 'reduced_library.fa.clstr.curated.xlsx' 

source ~/.bashrc

export PATH=$PATH:/data/alicel/chapter2/4_tools/seqtk

basePATH="/data/alicel/chapter2"
workdir="${basePATH}/1_dataTreatment/4_annotations_EarlGrey"

cd ${workdir}/2_single_library

# remove the ">" from the list (otherwise gives blank file)
sed -i 's/>//g' ${workdir}/2_single_library/REP_to_check.list
sed -i 's/>//g' ${workdir}/2_single_library/nonREP_to_check.list

# create a subset fasta file from the FASTA and the list
seqtk subseq ${workdir}/2_single_library/reduced_library.fa \
  ${workdir}/2_single_library/REP_to_check.list \
  > ${workdir}/2_single_library/subset_REP_seq.fa

seqtk subseq ${workdir}/2_single_library/complete_library.fa \
  ${workdir}/2_single_library/nonREP_to_check.list \
  > ${workdir}/2_single_library/subset_nonREP_seq.fa 