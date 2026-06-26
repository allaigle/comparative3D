#!/bin/bash

#SBATCH -c 1 # processor' number
#SBATCH --mem 16MB # memory

# [USAGE] sbatch -J D3 -c 1 --mem=16MB D3_rename_libraries.sh

basePATH="/data/alicel/chapter2"
workdir="${basePATH}/1_dataTreatment/4_annotations_EarlGrey"
listLibraries="/data/alicel/chapter2/0_data/0_collection/list_TE_libraries_species.txt" 
# 3cols: toby_lib species_lib species

mkdir -p ${workdir}/1_libraries_TBaril/renamed_libraries \
  ${workdir}/2_single_library

# cp and rename files
while read TOBY_LIB SPECIES_LIB _; do

  cp ${workdir}/1_libraries_TBaril/${TOBY_LIB} \
    ${workdir}/1_libraries_TBaril/renamed_libraries/${SPECIES_LIB}
    
done < $listLibraries 