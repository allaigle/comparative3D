#!/bin/bash
#
####################################################################################################################
#                              Correct raw.h5 matrices and normalize                         
#                        using HiCExplorer, bash and the SLURM scheduler                                                                             
#
# [Goal] Correct and normalize raw.h5 matrices using ICE normalization.
#
# [USAGE] sbatch -J E2 -c 1 --mem=16GB E2_correct_ICE_HiCExplorer.sh
#
# [Inputs] TXT file with thresholds & .raw.h5
# [Outputs] .corrected_ICE.h5
#
# [Adaptation & History]
# Mai-Aug 2025 Alice Laigle (alice.laigle@unine.ch) 
# Aug 2025 - AL with the help of Ashwini V. Mohan - reduced it to ICE, no KR
#
# [Documentation] 
# Correct: https://hicexplorer.readthedocs.io/en/latest/content/tools/hicCorrectMatrix.html
# Normalize: https://hicexplorer.readthedocs.io/en/latest/content/tools/hicNormalize.html
#
# NOTES: 
# 1. Thresholds have been determined using 'E1_hicpro2h5_diagnosticPlot_HiCExplorer.sh' script
# 3. Thresholds have been manually stored in a TXT containing 2 cols (\t): species xMinThreshold here: 
## '${basePATH}/0_data/0_collection/allDarwin_species_${shortRES}kbDiagThreshold.filtered_55.txt'
### where ${shortRES} can be either 10 or 50kb in our case - used only 50kb at the end.
###################################################################################################################
#
#SBATCH -c 1 # processor' number
#SBATCH --mem=16G

basePATH="/data/alicel/chapter2"
workdir="${basePATH}/1_dataTreatment/3_conversions/1_h5"
RESOLUTION="10000 50000" ## MIGHT BE CHANGED - would require the appropriate TXT file

mkdir -p ${workdir}/2_corrected_ICE

cd ${workdir}

conda activate hicexplorer

###################################################################################################################
### Correct and normalize matrices

for RES in $RESOLUTION; do 

   echo "Resolution: ${RES}"

  # Get the TXT file 
  shortRES=$(echo $(( RES / 1000 ))) 

  ## MIGHT BE CHANGED
  ## 2 cols (\t): species xMinThreshold
  INFILE="${basePATH}/0_data/0_collection/allDarwin_species_${shortRES}kbDiagThreshold.filtered_55.txt"

  echo "TXT file used as input: '${basePATH}/0_data/0_collection/allDarwin_species_${shortRES}kbDiagThreshold.filtered_55.txt'"

  ##### CORRECT

  while read SPECIES THRESHOLD; do

    echo ""
    echo "Correcting matrix for ${SPECIES} using '--correctionMethod ICE --filterThreshold ${THRESHOLD}' options for ${RES} resolution:"

    # ICE
    hicCorrectMatrix correct \
        --correctionMethod ICE \
        -m ${workdir}/1_raw/${SPECIES}_${RES}.raw.h5 \
        --filterThreshold ${THRESHOLD} 5 \
        -o ${workdir}/2_corrected_ICE/${SPECIES}_${RES}.corrected_ICE.h5   
  
    echo "Created: ${workdir}/2_corrected_ICE/${SPECIES}_${RES}.corrected_ICE.h5"

  done < $INFILE
done 

###################################################################################################################

conda deactivate

chmod -R 775 ${workdir} 

echo "Correction and normalization done."