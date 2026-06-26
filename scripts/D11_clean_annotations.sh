#!/bin/bash
#
####################################################################################################################
#                 Clean gene/TE annotations, convert and 
#          using BEDOPS, bedtools, bash and the SLURM scheduler                                                                             
#
# [GAOL] Clean annotations, convert to BED, filter gene annotation if overlaps with TE annotation and create
#        windows for each annotations to be plotted.
#
# [USAGE] sbatch -J D11 -c 1 --mem=8GB D11_clean_annotations.sh
#
# [Input] .fai, _braker.gff3 & .filteredRepeats.gff
# [Output] .gene_cov.bed and .TE_cov.bed
#
# [Adaptation & History]
# Nov 2024 Alice Laigle (alice.laigle@unine.ch)
# Feb 2025 AL - using BEDOPS to convert GFF (1-based) to BED (0-based)
# Mar 2025 AL - added TE cleaning steps & TE classes
#
# NOTE: BEDOPS and bedtools already installed.
###################################################################################################################

#SBATCH -c 1 # processor' number
#SBATCH --mem 8GB # memory

# Variables - TO BE CHANGED
source ~/.bashrc
basePATH="/data/alicel/chapter2"
workdir="${basePATH}/2_analysis/1_annotations"
INFILE="${basePATH}/0_data/0_collection/allDarwin_phylum_species_assembly_SRA.filtered_60.txt"

cd ${workdir}

# Create subdirectories
mkdir -p ${workdir}/0_inputs/gff3_genes \
  ${workdir}/0_inputs/cleaned_gff3_genes \
  ${workdir}/0_inputs/cleaned_gff_TEs \
  ${workdir}/0_inputs/fasta \
  ${workdir}/1_cleaned_genes \
  ${workdir}/2_BED_bedops \
  ${workdir}/3_subtract_anno/1_to_remove \
  ${workdir}/3_subtract_anno/2_filtered \
  ${workdir}/3_subtract_anno/3_summary \
  ${workdir}/4_make_windows_anno \
  ${workdir}/5_coverage_anno 


# INPUTS

while read -r PHYLUM SPECIES ASSEMBLY _; do
  
  # get 
  ln -sf ${basePATH}/0_data/1_genomes/wo500/${SPECIES}_${ASSEMBLY}.wo500.fasta.fai \
    ${workdir}/0_inputs/fasta/${SPECIES}.wo500.fasta.fai

  # get BRAKER gene anno
  ln -sf ${basePATH}/1_dataTreatment/2_annotations_braker/${SPECIES}/braker.gff3 \
    ${workdir}/0_inputs/gff3_genes/${SPECIES}_braker.gff3
  
  # get RepeatMasker TE anno
  ln -sf ${basePATH}/1_dataTreatment/2c_annotations_EarlGrey/3_final_annotation/${SPECIES}_EarlGrey/${SPECIES}_summaryFiles/${SPECIES}.filteredRepeats.gff \
    ${workdir}/0_inputs/cleaned_gff_TEs/${SPECIES}.filteredRepeats.gff

done < $INFILE	


#########################################################################################################################
# 1. Clean annotations

## 1.1. Genes
### 1.1.1. Keep only the chromosome name 
while read _ SPECIES _ _; do

  awk -F'\t' '{split($1, a, "_"); $1=a[1]}1' OFS='\t' ${workdir}/0_inputs/gff3_genes/${SPECIES}_braker.gff3 \
  > ${workdir}/0_inputs/cleaned_gff3_genes/${SPECIES}_braker.cleaned.gff3

done < ${INFILE}

### 1.1.2. Get only genes from GFF3 files
find ${workdir} -name "*_braker.cleaned.gff3*" | sort > ${workdir}/0_inputs/list_cleaned_gff3.txt # get full path & sort
list_cleaned_gff3="${workdir}/0_inputs/list_cleaned_gff3.txt"

while read FILE; do  
  SPECIES=$(basename "$FILE" | cut -d'_' -f1)

  grep "gene" $FILE > ${workdir}/1_cleaned_genes/${SPECIES}_braker.genes.gff3 

done < $list_cleaned_gff3 


## 1.2. TEs
### 1.2.1. Delete "ARTEFACT"

while read _ SPECIES _ _; do

  grep -v "ARTEFACT" ${workdir}/0_inputs/cleaned_gff_TEs/${SPECIES}.filteredRepeats.gff \
    > ${workdir}/0_inputs/cleaned_gff_TEs/${SPECIES}_filteredRepeats.cleaned.gff

  rm ${workdir}/0_inputs/cleaned_gff_TEs/${SPECIES}.filteredRepeats.gff

done < ${INFILE}


### 1.2.2. Merge as "Uncertain": "Low_complexity", "Simple_repeat", "SINE?"

PATTERNS="Low_complexity Simple_repeat SINE?"

while read _ SPECIES _ _; do

  for pattern in $PATTERNS; do 

    sed -i 's/'${pattern}'/Uncertain/g' ${workdir}/0_inputs/cleaned_gff_TEs/${SPECIES}_filteredRepeats.cleaned.gff
  
  done

done < ${INFILE}


#########################################################################################################################
# 2. Convert GFF to BED using BEDOPS

echo "Converting GFF to BED files for both gene and TE annotations."

while read -r _ SPECIES _ _; do

  # get chromosome name, start, end, "gene" or "TE" 
  gff2bed < ${workdir}/1_cleaned_genes/${SPECIES}_braker.genes.gff3 \
    > ${workdir}/2_BED_bedops/${SPECIES}.genes.sorted.bed

  gff2bed < ${workdir}/0_inputs/cleaned_gff_TEs/${SPECIES}_filteredRepeats.cleaned.gff \
    > ${workdir}/2_BED_bedops/${SPECIES}.TE.sorted.bed

done < $INFILE


#########################################################################################################################
# 3. Filter gene annotation if overlap more than 100bp with TEs

echo "Filtering gene annotation if it overlaps TE annotations for >100bp."

while read -r _ SPECIES _ _; do
  
  # Step 1: Get the list of gene lines to remove
  bedtools intersect -a ${workdir}/2_BED_bedops/${SPECIES}.genes.sorted.bed \
    -b ${workdir}/2_BED_bedops/${SPECIES}.TE.sorted.bed -wo | awk '$NF > 100' | cut -f1-3 \
    > ${workdir}/3_subtract_anno/1_to_remove/${SPECIES}.to_remove.bed

  # Step 2: Create a BED file (already sorted) of genes without those overlapping by 100bp with TEs
  bedtools subtract -a ${workdir}/2_BED_bedops/${SPECIES}.genes.sorted.bed \
    -b ${workdir}/3_subtract_anno/1_to_remove/${SPECIES}.to_remove.bed \
    > ${workdir}/3_subtract_anno/2_filtered/${SPECIES}.filtered_genes.bed

done < $INFILE


#############################################################################################################
# 4. Create 10kb window files of the genome

echo "Creating the <SPECIES>.wo500.windows10kb.bed files."

while read -r _ SPECIES _ _; do

  bedtools makewindows -g ${workdir}/0_inputs/fasta/${SPECIES}.wo500.fasta.fai -w 10000 \
    > ${workdir}/4_make_windows_anno/${SPECIES}.wo500.windows10kb.bed 

done < $INFILE	


#############################################################################################################
# 5. Get the gene/TE coverage per window

echo "Creating the <SPECIES>.gene_cov.bed and <SPECIES>.TE_cov.bed files."

while read -r _ SPECIES _ _; do

  # gene
  bedtools coverage -a ${workdir}/4_make_windows_anno/${SPECIES}.wo500.windows10kb.bed \
    -b ${workdir}/3_subtract_anno/2_filtered/${SPECIES}.filtered_genes.bed \
    > ${workdir}/5_coverage_anno/${SPECIES}.gene_cov.bed 

  # TE
  bedtools coverage -a ${workdir}/4_make_windows_anno/${SPECIES}.wo500.windows10kb.bed \
    -b ${workdir}/2_BED_bedops/${SPECIES}.TE.sorted.bed \
    > ${workdir}/5_coverage_anno/${SPECIES}.TE_cov.bed 
done < $INFILE

#############################################################################################################
chmod -R 775 ${workdir}