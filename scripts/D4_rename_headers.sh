#!/bin/bash

#SBATCH -c 1 # processor' number
#SBATCH --mem 16MB # memory

# [USAGE] sbatch -J D4 -c 1 --mem=16MB D4_rename_headers.sh

basePATH="/data/alicel/chapter2"
workdir="${basePATH}/1_dataTreatment/4_annotations_EarlGrey"

cd ${workdir}/1_libraries_TBaril/renamed_libraries

for file in *-families.fa.strained ;
	do species=${file%-*}
      sed -i 's/^>/>'${species}'___/g' ${file}
done
# note: "${species}'___" allowed me to clean it very easily later.