#!/bin/bash
#
####################################################################################################################
#                   Convert H5 matrices to ginteractions and create summary                    
#                       using HiCExplorer, bash and the SLURM scheduler                                                                             
#
# [USAGE] sbatch -J E3 -c 1 --mem=16GB E3_convert_h5_to_ginteractions_summarize_HiCExplorer.sh
#
# [Inputs] *_${RES}.${TYPE}.h5
# [Outputs] *_${RES}.${TYPE}.tsv
#
# [Adaptation & History]
# Aug 2025 Alice Laigle (alice.laigle@unine.ch) 
###################################################################################################################
#
#SBATCH -c=1
#SBATCH --mem=64MB

basePATH="/data/alicel/chapter2"
workdir="${basePATH}/1_dataTreatment/5_conversions"
INFILE="${basePATH}/0_data/0_collection/allDarwin_phylum_species_assembly_SRA.filtered_55.txt"
RESOLUTION="10000 50000" 
TYPE_DATA="raw corrected_ICE" 

cd ${workdir}

mkdir -p ${workdir}/2_ginteractions/1_raw \
  ${workdir}/2_ginteractions/2_corrected_ICE

###################################################################################################################
# 1. Convert H5 to GINTERACTIONS

conda activate hicexplorer

while read _ SPECIES _ _ ; do

  echo $SPECIES

  for RES in $RESOLUTION; do 

  hicConvertFormat \
    -m ${workdir}/1_h5/1_raw/${SPECIES}_${RES}.raw.h5 \
    --inputFormat h5 \
    --outputFormat ginteractions \
    -o ${workdir}/2_ginteractions/1_raw/${SPECIES}_${RES}.raw

  echo "ginteractions' conversion done for ${SPECIES} (raw) for ${RES} resolution."
  echo ""

  hicConvertFormat \
    -m ${workdir}/1_h5/2_corrected_ICE/${SPECIES}_${RES}${TYPE}h5 \
    --inputFormat h5 \
    --outputFormat ginteractions \
    -o ${workdir}/2_ginteractions/2_corrected_ICE/${SPECIES}_${RES}.corrected_ICE

  echo "ginteractions' conversion done for ${SPECIES} (corrected_ICE) for ${RES} resolution."
  echo ""

  done

done < $INFILE

###################################################################################################################
# 2. Create summary of pairs for raw and corrected_ICE matrices

for RES in $RESOLUTION; do

  echo $RES

  for TYPE in $TYPE_DATA; do

    echo $TYPE

    touch ${workdir}/allDarwin.summary_interactions.${RES}.${TYPE}.txt
    OUTPUTFILE="${workdir}/allDarwin.summary_interactions.${RES}.${TYPE}.txt"
    > $OUTPUTFILE # clean in case 

    # create header (8 columns)
    echo -e "Species\ttotal_interactions\tcis_interactions\tcis_perc\tshort_cis_interactions\tshort_cis_perc\tlong_cis_interactions\tlong_cis_perc\ttrans_interactions\ttrans_perc" \
      >> ${OUTPUTFILE}

    while read -r _ SPECIES _ _; do

      echo $SPECIES

      # get total number of interactions
      total_interactions=$(awk '{sum += $NF} END {printf "%.0f\n", sum}' ${workdir}/2_ginteractions/*/${SPECIES}_${RES}.${TYPE}.tsv)

      # count interactions
      cis_interactions=$(awk -F'\t' '$1 == $4 {sum += $NF} END {printf "%.0f\n", sum}' ${workdir}/2_ginteractions/*/${SPECIES}_${RES}.${TYPE}.tsv)
      trans_interactions=$(awk -F'\t' '$1 != $4 {sum += $NF} END {printf "%.0f\n", sum}' ${workdir}/2_ginteractions/*/${SPECIES}_${RES}.${TYPE}.tsv)
      
      # count short and long cis interactions
      short_cis_interactions=$(awk -F'\t' '$1 == $4 && ($2 == $5 || $2 == $5 + 50000) {sum += $NF} END {printf "%.0f\n", sum}' ${workdir}/2_ginteractions/*/${SPECIES}_${RES}.${TYPE}.tsv)
      long_cis_interactions=$(echo "$cis_interactions - $short_cis_interactions" | bc)

      # calculate percentages
      cis_perc=$(echo "scale=3; ($cis_interactions / $total_interactions) * 100" | bc)
      short_cis_perc=$(echo "scale=3; ($short_cis_interactions / $total_interactions) * 100" | bc)
      long_cis_perc=$(echo "scale=3; $cis_perc - $short_cis_perc" | bc)
      trans_perc=$(echo "scale=3; 100 - $cis_perc" | bc)
 
      # Print results
      printf "Total interactions: %.0f\n" "$total_interactions"
      printf "Cis interactions: %.0f (%.4f%%)\n" "$cis_interactions" "$cis_perc"
      printf "Cis short interactions: %.0f (%.4f%%)\n" "$short_cis_interactions" "$short_cis_perc"
      printf "Cis long interactions: %.0f (%.4f%%)\n" "$long_cis_interactions" "$long_cis_perc"
      printf "Trans interactions: %.0f (%.4f%%)\n" "$trans_interactions" "$trans_perc"

      ## Output result in the global tsv file
      echo -e "${SPECIES}\t${total_interactions}\t${cis_interactions}\t${cis_perc}\t${short_cis_interactions}\t${short_cis_perc}\t${long_cis_interactions}\t${long_cis_perc}\t${trans_interactions}\t${trans_perc}" \
        >> ${OUTPUTFILE}

    done < $INFILE

  done 

done 
