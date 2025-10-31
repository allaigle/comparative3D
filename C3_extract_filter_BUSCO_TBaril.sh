#!/bin/bash
#SBATCH -c 1 # processor' number
#SBATCH --mem 16MB # memory

# Adaptation Alice Laigle from Toby Baril's script.
# [USAGE] sbatch -J C3 -c 1 --mem=16MB C3_extract_filter_BUSCO_TBaril.sh

basePATH="/data/alicel/chapter2"
workdir="${basePATH}/1_dataTreatment/2_phyloTree"
INFILE="${basePATH}/0_data/0_collection/phylotree_phylum_species_assembly.63.txt"

## Move to parent directory as this is where we will store the list
cd ${workdir}

## 1. Extract any gene from each directory that is complete
for file in ${workdir}/1_busco/*/run_fungi_odb10/full_table.tsv ;
	do grep -v "^#" ${file} | awk '$2 == "Complete" {print $1}' >> complete_busco_ids.txt ;
done

for file in ${workdir}/1_busco/*/run_fungi_odb10/full_table.tsv ;
	do echo $file
done

var=$(wc -l complete_busco_ids.txt | cut -d' ' -f1) # 45485 
echo "Number of complete_busco_ids.txt: $var"

## 2. Sort the complete IDs and count how many genomes they appeared in
sort complete_busco_ids.txt | uniq -c > complete_busco_ids_with_counts.txt

## 3. Keep only those present in at least 3 genomes
awk '$1 > 2 {print $2}' complete_busco_ids_with_counts.txt > final_busco_ids.txt
var=$(wc -l final_busco_ids.txt | cut -d' ' -f1) # 758 
echo "Number of final_busco_ids.txt: $var"