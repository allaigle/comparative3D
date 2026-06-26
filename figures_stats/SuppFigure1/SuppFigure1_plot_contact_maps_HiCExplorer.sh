#!/bin/bash
#
####################################################################################################################
#                               Plot corrected iced matrices                         
#                        using HiCExplorer, bash and the SLURM scheduler                                                                             
#
# [Goal] Plot raw and corrected_ICE matrices from HiCExplorer for whole genome and per chromosome (only corrected).
#
# [USAGE] sbatch -J plot_maps -c 1 --mem=16GB SuppFigure1_plot_contact_maps_HiCExplorer.sh
#
# [Input] .corrected_ICE.h5
# [Output] .png for corrected_ICE at different resolutions and vMax
#
# [Adaptation & History]
# Aug 2025 Alice Laigle (alice.laigle@unine.ch) with the help of Ashwini V. Mohan
#
# [Documentation] https://hicexplorer.readthedocs.io/en/latest/content/tools/hicPlotMatrix.html
#
# [Notes]
# 1. HiCExplorer is in an environement that has to be activated before running the sbatch command
# 2. Tested different vMax for each resolutions to see which vMax is the most appropriate for this dataset.
###################################################################################################################
#
#SBATCH -c=1
#SBATCH --mem=16G

# Variables
basePATH="/data/alicel/chapter2"
workdir="${basePATH}/1_dataTreatment/5_conversions/1_h5"
INFILE="${basePATH}/0_data/0_collection/allDarwin_phylum_species_assembly_SRA.filtered_55.txt"
RESOLUTION="10000 50000" # MIGHT BE CHANGED - used 5000 at the end
vMax="1000 2500 5000" # MIGHT BE CHANGED - used 2500 at the end

mkdir -p ${workdir}/4_contact_maps

cd ${workdir}

conda activate hicexplorer

###################################################################################################################
# PLOT
for RES in $RESOLUTION; do 

  echo "Resolution: ${RES}"

  while read _ SPECIES _ _; do

    echo "Species: ${SPECIES}"
    mkdir -p ${workdir}/4_contact_maps/${RES}

    for VMAX in $vMax; do

      echo $VMAX

      #### Corrected ICE
      echo ""
      echo "Plotting corrected_ICE matrix for ${SPECIES} using '--log1p --vMax ${VMAX}' options for ${RES} resolution:"

      hicPlotMatrix \
          --matrix ${workdir}/2_corrected_ICE/${SPECIES}_${RES}.corrected_ICE.h5 \
          --outFileName ${workdir}/4_contact_maps/${RES}/${SPECIES}_${RES}.corrected_ICE.log1p_vMax${VMAX}.png \
          --log1p \
          --vMax $VMAX \
          --rotationX 90 \
          --dpi 300

      echo ""
      echo "Created: ${workdir}/4_contact_maps/${RES}/${SPECIES}_${RES}.corrected_ICE.log1p_vMax${VMAX}.png"

    done 
  done < $INFILE    
done 

conda deactivate

chmod -R 775 ${workdir} 

echo "Plotting contact maps done."