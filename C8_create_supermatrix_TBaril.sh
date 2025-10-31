#!/bin/bash

#SBATCH -c 1 # processor' number
#SBATCH --mem 32MB # memory

# Adaptation Alice Laigle from Toby Baril's script.
# [USAGE] sbatch -J C8 -c 1 --mem 32MB C8_create_supermatrix_TBaril.sh

basePATH="/data/alicel/chapter2"
workdir="/data/alicel/chapter2/1_dataTreatment/2_phyloTree"

## make a directory
mkdir -p ${workdir}/4_supermatrix && cd ${workdir}/4_supermatrix

## make a mamba env with the correct tools - conda in my case (miniconda3 to be precise)
#conda create -n phykit bioconda::phykit bioconda::clipkit bioconda::iqtree 
## - got error, due to biopython version, gonna do step by step
#conda create -n phykit
conda activate phykit
#conda install bioconda
#conda install bioconda::phykit bioconda::clipkit bioconda::iqtree

## create a list of alignment files
realpath ${workdir}/3_alignments/fixedHeader/* > alignmentList.txt

## concatenate the alignments
phykit create_concat -a alignmentList.txt -p fungi_supermatrix
## made a note of the stout as phykit_log.txt

## trim uninformative regions
clipkit fungi_supermatrix.fa
## made a note of the stout as clipkit_log.txt

chmod -R 775 ${workdir}/4_supermatrix