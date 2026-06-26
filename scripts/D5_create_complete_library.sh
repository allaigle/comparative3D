#!/bin/bash

#SBATCH -c 1 # processor' number
#SBATCH --mem 32MB # memory

# [USAGE] sbatch -J D5 -c 1 --mem=32MB D5_create_complete_library.sh

basePATH="/data/alicel/chapter2"
workdir="${basePATH}/1_dataTreatment/4_annotations_EarlGrey"

for file in ${workdir}/1_libraries_TBaril/renamed_libraries/*-families.fa.strained ; 
  do 
    cat $file >> ${workdir}/2_single_library/complete_library.fa
done