#!/bin/bash

#SBATCH -c 1 # processor' number
#SBATCH --mem 2MB # memory

# Adaptation Alice Laigle from Toby Baril's script.
# [USAGE] sbatch -J C9 -c 1 --mem 2MB C9_run_iqtree2_TBaril.sh

basePATH="/data/alicel/chapter2"
workdir="/data/alicel/chapter2/1_dataTreatment/2_phyloTree"

## make a directory for the tree
mkdir -p ${workdir}/5_tree && cd ${workdir}/5_tree

## symlink the alignment that the tree will use
ln -sr ${workdir}/4_supermatrix/fungi_supermatrix.fa.clipkit .

# Recommended by Toby:
sbatch -c 16 -p normal.1000h --mem=64G \
  --wrap="iqtree2 -s fungi_supermatrix.fa.clipkit -bb 1000 -T 16 \
  -o Cmojavensis,Cobscurus,Lpennispora --verbose"

#NOTES: 
# -s : Specify input alignment file in PHYLIP, FASTA, NEXUS, CLUSTAL or MSF format.
# -bb : Specify number of bootstrap replicates (>=1000).
# -T : Specify the number of CPU cores to use only for the SH-aLRT test. 
## If -T AUTO is specified, IQ-TREE will use all available cores. 
## NOTE: this option has no effect on tree search, which is still single-threaded.
# -o : Specify an outgroup taxon name to root the tree. The output tree in .treefile will be rooted accordingly. 
## DEFAULT: first taxon in alignment
