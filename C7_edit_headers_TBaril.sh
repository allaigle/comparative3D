#!/bin/bash

#SBATCH -c 1 # processor' number
#SBATCH --mem 4MB # memory

# Adaptation Alice Laigle from Toby Baril's script.
# [USAGE] sbatch -J C7 -c 1 --mem=4MB C7_edit_headers_TBaril.sh

basePATH="/data/alicel/chapter2"
workdir="/data/alicel/chapter2/1_dataTreatment/2_phyloTree"

cd $workdir 

## make a directory to store results (rather than overwrite in case we mess something up!)
mkdir -p ${workdir}/3_alignments/fixedHeader 

## delete empty alignments
find . -size  0 -print -delete

## MAKE A NOTE OF WHICH ARE PRINTED!
## - none of them was empty

## fix the fasta headers
for file in *.aln ;
	do sed '/^>/s/|.*//g' $file > ${workdir}/3_alignments/fixedHeader/${file} ; 
done