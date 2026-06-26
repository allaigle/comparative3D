#!/bin/bash
#
####################################################################################################################
#                  Convert raw Hi-C matrix to h5 and run diagnostic plot              
#                    using HiCExplorer, bash and the SLURM scheduler                                                                             
#
# [USAGE] sbatch -J E1 -c 1 --mem=16GB E1_hicpro2h5_diagnosticPlot_HiCExplorer.sh
#
# [Input] allDarwin_phylum_species_assembly_SRA.filtered_55.txt, .matrix & _abs.bed
# [Output] raw.h5 & .diag_raw.pdf
#
# [Adaptation & History]
# Oct-Nov 2024 Alice Laigle (alice.laigle@unine.ch)
# April 2025 AL - added Fit-Hi-C & ginteractions adaptations
# May 2025 AL - added merged_nodups short format with score 
# May 2025 AL - splitted script as added diagnostic_plot and modified normalization steps
# Aug 2025 AL - simplified it to .h5, 50 kb and diag plots (rest unnecessary now or after correction)
#
# NOTES:
#    1. 3DGB is already done
#    2. HiCExplorer are already installed, HiCExplorer having his own 'hicexplorer' conda env
#
# HiCExplorer documentation: https://hicexplorer.readthedocs.io/en/latest/
###################################################################################################################
#
#SBATCH -J E1
#SBATCH -c 1 # processor' number
#SBATCH --mem 16GB # memory

###################################################################################################################
#### ENVIRONMENT 

basePATH="/data/alicel/chapter2"
workdir="${basePATH}/1_dataTreatment/3_conversions"
RES="50000" # also tested 10, 20, 30, 40 & 100 kb

## 4 cols (\t): phylum species assembly SRA 
### Only necessary: phylum species SRA 
INFILE="${basePATH}/0_data/0_collection/allDarwin_phylum_species_assembly_SRA.filtered_55.txt"

mkdir -p ${workdir}/0_inputs/raw_abs \
  ${workdir}/0_inputs/raw_matrix \
  ${workdir}/1_h5/0_diagnostic_plots \
  ${workdir}/1_h5/1_raw

cd ${workdir}

###################################################################################################################
#### INPUTS

while read PHYLUM SPECIES _ SRA; do

  ln -sf ${basePATH}/1_dataTreatment/1_3DGB/species/${PHYLUM}/${SPECIES}/HiC-Pro/matrix/${SRA}/raw/${RES}/${SRA}_${RES}.matrix \
    ${workdir}/0_inputs/raw_matrix/${SPECIES}_${RES}.matrix

  ln -sf ${basePATH}/1_dataTreatment/1_3DGB/species/${PHYLUM}/${SPECIES}/HiC-Pro/matrix/${SRA}/raw/${RES}/${SRA}_${RES}_abs.bed \
    ${workdir}/0_inputs/raw_abs/${SPECIES}_${RES}_abs.bed

done < $INFILE

##################################################################################################################
# 1. Get .h5 matrix from raw data from HiC-Pro

cd ${workdir}

conda activate hicexplorer

while read _ SPECIES _ _; do

  echo $SPECIES
  
  hicConvertFormat -m ${workdir}/0_inputs/raw_matrix/${SPECIES}_${RES}.matrix \
    --bedFileHicpro ${workdir}/0_inputs/raw_abs/${SPECIES}_${RES}_abs.bed \
    --inputFormat hicpro \
    --outputFormat h5 \
    -o ${workdir}/1_h5/1_raw/${SPECIES}_${RES}.raw.h5

  echo ""
  echo "H5 conversions done for ${SPECIES}."

done < $INFILE

##################################################################################################################
# 2. Diagnostic plot for raw.h5 matrix

while read _ SPECIES _ _; do

  echo $SPECIES

  hicCorrectMatrix diagnostic_plot \
    --matrix ${workdir}/1_h5/1_raw/${SPECIES}_${RES}.raw.h5 \
    -o ${workdir}/1_h5/0_diagnostic_plots/${SPECIES}_${RES}.diag_raw.pdf \
    --verbose

  echo ""
  echo "Diagnostic plots done for ${SPECIES}."

done < $INFILE

##################################################################################################################
# Cleaning

conda deactivate

chmod -R 775 ${workdir}