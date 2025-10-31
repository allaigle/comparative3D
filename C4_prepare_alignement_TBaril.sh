#!/bin/bash
#SBATCH -c 1 # processor' number
#SBATCH --mem 16MB # memory

# Adaptation Alice Laigle from Toby Baril's script.
# [USAGE] sbatch -J C4 -c 1 --mem=16MB C4_prepare_alignement_TBaril.sh

basePATH="/data/alicel/chapter2"
workdir="${basePATH}/1_dataTreatment/2_phyloTree"
INFILE="${basePATH}/0_data/0_collection/phylotree_phylum_species_assembly.63.txt"

cd ${workdir}

# 1. Copy amino acid sequences and name files with species name in them

## make a directory to store amino sequences 
mkdir -p ${workdir}/2_busco_aa

## extract and rename amino acid sequences
# Note: had to make some modififactions to extract the species name from Toby's notes
#sedCommand1="s|/data/alicel/chapter2/1_dataTreatment/2_phyloTree/1_busco/BUSCO_||g"
#sedCommand2="s|/run_fungi_odb10/busco_sequences/single_copy_busco_sequences||g"

for DIRECTORY in $(find ${workdir}/1_busco/ -type d -name "single_copy_busco_sequences") ; 
	do SPECIES=$(echo ${DIRECTORY} | sed -e 's|/data/alicel/chapter2/1_dataTreatment/2_phyloTree/1_busco/BUSCO_||g;s|/run_fungi_odb10/busco_sequences/single_copy_busco_sequences||g' ) ;
		for FILE in ${DIRECTORY}/*.faa ;
			do filename=$(basename ${FILE}) ; 
        echo $SPECIES
        echo $filename
        cp ${FILE} ${workdir}/2_busco_aa/${SPECIES}_${filename} ;
		done ;
done

# 2. Add species names to FASTA headers

## move into the correct directory
cd ${workdir}/2_busco_aa

## Add names to files
for file in *.faa ; 
	do sed -i 's/^>/>'${file%_*}'|/g; /^>/s/:.*//g' ${file} ;
done

# 3. Generate a single file for each BUSCO ID

while read line ;
	do cat *_${line}.faa >> ${line}_aa.fasta ;
done < ${workdir}/final_busco_ids.txt

chmod -R 775 ${workdir}