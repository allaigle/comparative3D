#!/bin/bash

#SBATCH -c 1 # processor' number
#SBATCH --mem 16MB # memory

# Adaptation Alice Laigle from Toby Baril's script.
# [USAGE] sbatch -J D8 -c 1 --mem=16MB D8_count_N_nt.list_to_remove.sh

# for each seq, count the number of "N" nt at the beginning and at the end, 
# and for the whole sequence
## > If 20% of the 100nt (50 first & 50 last)
## > If 10% of whole sequence
## >> Delete from the reduced and pick the second choice

source ~/.bashrc
export PATH=$PATH:/data/alicel/chapter2/4_tools/seqtk
basePATH="/data/alicel/chapter2"
workdir="${basePATH}/1_dataTreatment/4_annotations_EarlGrey"
cd ${workdir}/2_single_library

TYPE_SEQ="REP   nonREP" # REP or nonREP

for TYPE in $TYPE_SEQ; do

  # 1. get the first and last 50 nt of each seq - do same for nonREP
  # 1.1. 50 first
  seqkit subseq -r 1:50 ${workdir}/2_single_library/subset_${TYPE}_seq.fa \
    > ${workdir}/2_single_library/head50_subset_${TYPE}_seq.fa 

  # 1.2. 50 last
  seqkit subseq -r -50:-1 ${workdir}/2_single_library/subset_${TYPE}_seq.fa \
    > ${workdir}/2_single_library/tail50_subset_${TYPE}_seq.fa # 50 last; do same for REP 

  # 2. combine them
  seqkit concat ${workdir}/2_single_library/head50_subset_${TYPE}_seq.fa \
    ${workdir}/2_single_library/tail50_subset_${TYPE}_seq.fa \
    > ${workdir}/2_single_library/subset100nt_${TYPE}_seq.fa 

  # 3. get the nb of "N" for each sequence 
  seqtk comp ${workdir}/2_single_library/subset100nt_${TYPE}_seq.fa  \
    | awk -v OFS='\t' '{x=$3+$4+$5+$6;y=$2;print $1,y,y-x }' \
    > ${workdir}/2_single_library/count_N_nt_subset100nt_${TYPE}_seq.txt
  # result: Name of seq, nb nt of seq (100), nb N of seq (= % of N of seq)
  # e.g., Rirregularis___rnd-1_family-18#Unknown	100	2

  # 4. get only lines if >=20 of "N" to remove them
  ## >> if more than 20, remove the seq from REP
  ## & take the second possibility in the file
  awk '( $3 >= 20 )' ${workdir}/2_single_library/count_N_nt_subset100nt_${TYPE}_seq.txt \
    > ${workdir}/2_single_library/list_${TYPE}_seq_to_delete.100nt.txt 

done


# 5. Check if "N" represent more than 10% of the whole seq 
## Output contains 4 cols: seqName; nb nt; nb nt minus nb N; % of N nt in the seq  

for TYPE in $TYPE_SEQ; do 

  awkCommand='{x=$3+$4+$5+$6; y=$2; print $1, y, y-x, (y-x)/y*100}'

  seqtk comp ${workdir}/2_single_library/subset_${TYPE}_seq.fa \
   | awk -v OFS='\t' "$awkCommand" \
    > ${workdir}/2_single_library/count_N_nt_subset_${TYPE}_seq.txt

  awk '($4 >= 10)' ${workdir}/2_single_library/count_N_nt_subset_${TYPE}_seq.txt \
    > ${workdir}/2_single_library/list_${TYPE}_seq_to_delete.wholeSeq.txt 

done


# 6. Concatenate the list_${TYPE}_seq_to_delete.txt and 
## final_list_${TYPE}_seq_to_delete.txt and keep uniq seqNames

for TYPE in $TYPE_SEQ; do

  cat ${workdir}/2_single_library/list_${TYPE}_seq_to_delete.wholeSeq.txt \
    ${workdir}/2_single_library/list_${TYPE}_seq_to_delete.100nt.txt \
    > ${workdir}/2_single_library/tmp.${TYPE}.concat.txt

  awk '!seen[$1]++ { print $0 }' ${workdir}/2_single_library/tmp.${TYPE}.concat.txt \
  > ${workdir}/2_single_library/final_list_${TYPE}_seq_to_delete.txt

done

rm tmp*