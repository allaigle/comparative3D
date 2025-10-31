#!/bin/bash
#
####################################################################################################################
#                             Clean successful 3DGB workdir                         
#                            using bash and the SLURM scheduler                                                                             
#
# [USAGE] sbatch -c1 --mem=2MB B3_clean_3DGB_workdir.sh PHYLUM SPECIES SRA
#
# [Adaptation & History]
# Oct 2024 Alice Laigle (alice.laigle@gmail.com)
###################################################################################################################
#
#SBATCH -c 1 # processor' number
#SBATCH --mem 2MB # memory

PHYLUM=$1 # e.g., PHYLUM="basidio"
SPECIES=$2 # e.g., SPECIES="Abisporus"
SRA=$3 # e.g., "ERR9580487"

# Variable 
basePATH="/data/alicel/chapter2" # could be added as argument or modify, but not necessary in my case
workdir="${basePATH}/1_dataTreatment/1_3DGB"

cd ${workdir}

# Clean HiC-Pro subfolder 
mv ${SPECIES}/HiC-Pro/output/hic_results/data/${SRA}/${SRA}.allValidPairs ${SPECIES}/${SPECIES}_${SRA}.allValidPairs
mv ${SPECIES}/HiC-Pro/output/hic_results/pic ${SPECIES}/HiC-Pro/ ;\
mv ${SPECIES}/HiC-Pro/pic/${SRA}/* ${SPECIES}/HiC-Pro/pic/ ; rm -r ${SPECIES}/HiC-Pro/pic/${SRA} 
mv ${SPECIES}/HiC-Pro/output/hic_results/stats ${SPECIES}/HiC-Pro/ ;
mv ${SPECIES}/HiC-Pro/stats/${SRA}/* ${SPECIES}/HiC-Pro/stats/ ; rm -r ${SPECIES}/HiC-Pro/stats/${SRA} 
mv ${SPECIES}/HiC-Pro/output/logs ${SPECIES}/HiC-Pro
mv ${SPECIES}/HiC-Pro/output/config.txt ${SPECIES}/HiC-Pro
mv ${SPECIES}/HiC-Pro/output/hic_results/matrix ${SPECIES}/HiC-Pro

# Delete "useless" HiC-Pro subfolders  
rm -r ${SPECIES}/HiC-Pro/output
rm -r ${SPECIES}/HiC-Pro/bowtie2_index/
rm -r ${SPECIES}/HiC-Pro/merged_output
rm -r ${SPECIES}/HiC-Pro/merged_samples

# Delete other 3DGB generated folders and FASTQ files
rm -r ${SPECIES}/fastq_files
rm -r ${SPECIES}/sequence
rm -r ${SPECIES}/.snakemake
rm -r ${SPECIES}/logs

# Clean structure folder by keeping only cleaned ones
rm ${SPECIES}/structure/*/structure_completed.pdb
rm ${SPECIES}/structure/*/structure_verified_contigs.pdb
rm ${SPECIES}/structure/*/structure_with_chr.pdb

# Move the "SPECIES" folder to the finished SPECIES/PHYLUM one
mv ${SPECIES} species/${PHYLUM}
chmod -R 775 species/${PHYLUM}/${SPECIES} 
