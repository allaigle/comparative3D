#!/bin/bash

#SBATCH -c 1 # processor' number
#SBATCH --mem 4MB # memory

# Adaptation Alice Laigle from Toby Baril's script.
# [USAGE] sbatch -J C5 -c 1 --mem=4MB C5_prepare_run_mafft_TBaril.sh

basePATH="/data/alicel/chapter2"
workdir="/data/alicel/chapter2/1_dataTreatment/2_phyloTree"

## create mafft env and install it
#conda create -n mafft
conda activate mafft
#conda install conda-forge::mafft

## create list of the AA fasta files 
cd ${workdir}/2_busco_aa 
ls *_aa.fasta > ../alignmentParameters.txt 

## make a directory for the results
mkdir -p ${workdir}/3_alignments && cd ${workdir}/3_alignments

ln -s ${basePATH}/3_scripts/C6_run_mafft_TBaril.sh C6_run_mafft_TBaril.sh

## run mafft
bash C6_run_mafft_TBaril.sh

## get list of complete aligned IDs
ls *.aln > ../complete_aligned_ids.txt

## clean
mkdir -p ${workdir}/3_alignments/slurm ; mv slurm-* slurm/
chmod -R 775 ${workdir}/*