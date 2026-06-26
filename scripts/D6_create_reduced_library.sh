#!/bin/bash

#SBATCH -p normal.1000h # long queue
#SBATCH -c 1 # processor' number
#SBATCH --mem 16MB # memory

# Create consensus library

# [USAGE] sbatch -J D6 -c 1 -p normal.1000h --mem=16MB D6_create_reduced_library.sh

# Requires CD-HIT (v4.8.1; Li and Godzik, 2006)
## https://academic.oup.com/bioinformatics/article/22/13/1658/194225?login=false

source ~/.bashrc
basePATH="/data/alicel/chapter2"
workdir="${basePATH}/1_dataTreatment/4_annotations_EarlGrey"

cd ${workdir}/2_single_library

conda activate te_annot # created an env with some tools for TE annotations

cd-hit-est -i ${workdir}/2_single_library/complete_library.fa \
  -o ${workdir}/2_single_library/reduced_library.fa \
  -d 0 -aS 0.8 -c 0.8 -G 0 -g 1 -b 500 -T 32 -M 16000