#!/bin/bash

#SBATCH -c 1 # processor' number
#SBATCH --mem 16MB # memory

# [USAGE] sbatch -J D9 -c 1 --mem=16MB D9_create_final_curated_single_library.sh

# From the complete library, get all the sequences of the curated single library.

source ~/.bashrc
export PATH=$PATH:/data/alicel/chapter2/4_tools/seqtk
basePATH="/data/alicel/chapter2"
workdir="${basePATH}/1_dataTreatment/4_annotations_EarlGrey"

cd ${workdir}/2_single_library

sed -i 's/>//g' ${workdir}/2_single_library/final_list.reduced_library.fa.txt

seqtk subseq ${workdir}/2_single_library/complete_library.fa \
  ${workdir}/2_single_library/final_list.reduced_library.fa.txt \
  > ${workdir}/2_single_library/reduced_library.curated.fa

# Re-format the curated FASTA
seqtk subseq ${workdir}/2_single_library/complete_library.fa \
  ${workdir}/2_single_library/final_list.reduced_library.fa.txt \
  > ${workdir}/2_single_library/tmp.fa

# Format to 50 nt-long lines
seqkit seq -w 50 ${workdir}/2_single_library/tmp.fa \
  > ${workdir}/2_single_library/reduced_library.curated.fa

rm tmp.fa &&  chmod -R 775 *