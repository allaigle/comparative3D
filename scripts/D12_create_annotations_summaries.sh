#!/bin/bash
#
####################################################################################################################
#                      Create summary files 
#                 using bash and the SLURM scheduler                                                                             
#
# [GAOL] Create summaries for future plots.
#
# [USAGE] sbatch -J D12 -c 1 --mem=8MB D12_create_annotations_summaries.sh
#
# [Input] .fa, .fai, _chromSize.txt, _braker.gff3 & .filteredRepeats.gff
# [Output] .bed, .txt & .bedgraph of GC contents
#
# [Adaptation & History]
# Nov 2024 Alice Laigle (alice.laigle@unine.ch)
# Feb 2025 AL - adapted to new TE annotation
# Feb 2025 AL - using BEDOPS to convert GFF (1-based) to BED (0-based)
# Mar 2025 AL - added TE cleaning steps & TE classes
#
###################################################################################################################

#SBATCH -c 1 # processor' number
#SBATCH --mem 8GB # memory

# Variables - TO BE CHANGED
source ~/.bashrc
basePATH="/data/alicel/chapter2"
workdir="${basePATH}/2_analysis/1_annotations"
INFILE="${basePATH}/0_data/0_collection/allDarwin_phylum_species_assembly_SRA.filtered_60.txt"

cd ${workdir}

mkdir -p ${workdir}/0_inputs/chrSize \
  ${workdir}/6_TE_classes \
  ${workdir}/7_global_summaries


# INPUTS

while read -r PHYLUM SPECIES ASSEMBLY _; do
  # get fasta & its .fai
  ln -sf ${basePATH}/0_data/1_genomes/wo500/${SPECIES}_${ASSEMBLY}.fasta \
    ${workdir}/0_inputs/fasta/${SPECIES}.fasta

  ln -sf ${basePATH}/0_data/1_genomes/wo500/${SPECIES}_${ASSEMBLY}.fasta.fai \
    ${workdir}/0_inputs/fasta/${SPECIES}.fasta.fai

  # get chromSize 
  ln -sf ${basePATH}/1_dataTreatment/1_3DGB/species/${PHYLUM}/${SPECIES}/HiC-Pro/chromosome_sizes.txt \
    ${workdir}/0_inputs/chrSize/${SPECIES}_chromSize.txt

done < $INFILE	

#########################################################################################################################
# 1. Summary of gene annotation's filtering

touch ${workdir}/7_global_summaries/allSpecies.summary.gene_filtered.txt
OUTPUTFILE="${workdir}/7_global_summaries/allSpecies.summary.gene_filtered.txt"
> $OUTPUTFILE # clean in case 

# create header (8 columns)
echo -e "Species\tGenomeSize_bp_wo500\tGene_number\tGene_cov_perc\tGene_number_filtered\tGene_cov_filtered_perc\tDiff_number\tDiff_perc" \
  >> ${OUTPUTFILE}

while read -r _ SPECIES _ _; do

  ## get genome size after trimming
  Gsize_wo500=$(awk '{sum += $2} END {print sum}' ${workdir}/0_inputs/fasta/${SPECIES}.wo500.fasta.fai)
  
  ## calculate the number and the total length
  ### genes
  gene_nb=$(wc -l ${workdir}/2_BED_bedops/${SPECIES}.genes.sorted.bed | cut -d ' ' -f 1)
  gene_cov=$(awk '{sum += $3 - $2} END {print sum}' ${workdir}/2_BED_bedops/${SPECIES}.genes.sorted.bed)

  ### filtered genes
  filter_gene_nb=$(wc -l ${workdir}/3_subtract_anno/2_filtered/${SPECIES}.filtered_genes.bed | cut -d ' ' -f 1)
  filter_gene_cov=$(awk '{sum += $3 - $2} END {print sum}' ${workdir}/3_subtract_anno/2_filtered/${SPECIES}.filtered_genes.bed)

  ## calculate percentage

  gene_perc=$(echo "scale=3; ($gene_cov * 100) / $Gsize_wo500" | bc)
  filter_gene_perc=$(echo "scale=3; ($filter_gene_cov * 100) / $Gsize_wo500" | bc)

  ## calculate difference
  diff_gene_nb=$(($gene_nb-$filter_gene_nb))
  diff_gene_nb_perc=$(echo "${gene_perc}-${filter_gene_perc}" | bc)
  
  ## Output result in a summary file
  echo "Species: ${SPECIES}." \
    >> ${workdir}/3_subtract_anno/3_summary/${SPECIES}.summary.gene_filtered.txt
  echo "Total genome size after removing the <500kb scaffolds/contigs/chr: ${Gsize_wo500} bp." \
    >> ${workdir}/3_subtract_anno/3_summary/${SPECIES}.summary.gene_filtered.txt 
  echo "Number of annotated genes: ${gene_nb}." \
    >> ${workdir}/3_subtract_anno/3_summary/${SPECIES}.summary.gene_filtered.txt
  echo "Percentage of gene coverage: ${gene_perc}%." \
    >> ${workdir}/3_subtract_anno/3_summary/${SPECIES}.summary.gene_filtered.txt
  echo "Number of filtered genes: ${filter_gene_nb}." \
    >> ${workdir}/3_subtract_anno/3_summary/${SPECIES}.summary.gene_filtered.txt
  echo "Percentage of filtered gene coverage: ${filter_gene_perc}%." \
    >> ${workdir}/3_subtract_anno/3_summary/${SPECIES}.summary.gene_filtered.txt
  echo "Difference of gene number and filtered gene number: ${diff_gene_nb}." \
    >> ${workdir}/3_subtract_anno/3_summary/${SPECIES}.summary.gene_filtered.txt
  echo "Difference of gene cov percentage and filtered gene cov percentage: ${diff_gene_nb_perc}%." \
    >> ${workdir}/3_subtract_anno/3_summary/${SPECIES}.summary.gene_filtered.txt

  ## Output result in the global tsv file
  echo -e "${SPECIES}\t${Gsize_wo500}\t${gene_nb}\t${gene_perc}\t${filter_gene_nb}\t${filter_gene_perc}\t${diff_gene_nb}\t${diff_gene_nb_perc}" \
    >> ${OUTPUTFILE}

done < $INFILE


#############################################################################################################

## 2. Get the metrics

touch ${workdir}/7_global_summaries/allDarwin.summary.percent_Anno.txt
OUTPUTFILE="${workdir}/7_global_summaries/allDarwin.summary.percent_Anno.txt"
> ${workdir}/7_global_summaries/allDarwin.summary.percent_Anno.txt # clear in case

# create header (11 columns)
echo -e "Species\tGenomeSize_bp\tGenomeSize_bp_wo500\tChromosome_Number\tratioChrNbChrLength\tGene_Number\tGene_cov\tGene_cov_perc\tTE_Number\tTE_cov\tTE_cov_perc" \
  >> ${OUTPUTFILE}

while read -r _ SPECIES _ _; do
  
  ### 2.2.1. get genome size before and after trimming
  Gsize=$(awk '{sum += $2} END {print sum}' ${workdir}/0_inputs/fasta/${SPECIES}.fasta.fai)
  Gsize_wo500=$(awk '{sum += $2} END {print sum}' ${workdir}/0_inputs/fasta/${SPECIES}.wo500.fasta.fai)

  ### 2.2.2. get number chromosome & ratio of chromosome number/mean of chromosome length

  numberChr=$(wc -l <  ${workdir}/0_inputs/chrSize/${SPECIES}_chromSize.txt)

  meanChrLength=$(awk '{ total += $2 } END { print total/NR }'  ${workdir}/0_inputs/chrSize/${SPECIES}_chromSize.txt)

  ratioChrNbLength=$(awk "BEGIN {printf \"%.10f\", $numberChr / $meanChrLength}")

  ### 2.2.3. calculate the number and the total length
  # genes
  geneNumber=$(wc -l ${workdir}/3_subtract_anno/2_filtered/${SPECIES}.filtered_genes.bed| cut -d' ' -f1) 

  geneCov=$(awk '{sum += $3 - $2} END {print sum}' ${workdir}/3_subtract_anno/2_filtered/${SPECIES}.filtered_genes.bed)

  # TE
  TEnumber=$(wc -l ${workdir}/2_BED_bedops/${SPECIES}.TE.sorted.bed | cut -d' ' -f1) 

  TECov=$(awk '{sum += $3 - $2} END {print sum}' ${workdir}/2_BED_bedops/${SPECIES}.TE.sorted.bed)

  ### 2.2.4. Calculate percentage

  percentageGene=$(echo "scale=3; ($geneCov * 100) / $Gsize_wo500" | bc)
  percentageTE=$(echo "scale=3; ($TECov * 100) / $Gsize_wo500" | bc)

  # Output result in a summary file
  echo "Species: ${SPECIES}" >> ${workdir}/5_coverage_anno/${SPECIES}.percent_Anno.txt
  echo "Total Genome Size (bp): ${Gsize}" >> ${workdir}/5_coverage_anno/${SPECIES}.percent_Anno.txt
  echo "Total Genome Size after removing the <500kb scaffolds/contigs/chr (bp): ${Gsize_wo500}" \
    >> ${workdir}/5_coverage_anno/${SPECIES}.percent_Anno.txt
  echo "Chromosome number from Genome Size without <500kb: ${numberChr}" \
    >> ${workdir}/5_coverage_anno/${SPECIES}.percent_Anno.txt 
  echo "Ratio chromosome Number/mean chromosome length (Genome Size without <500kb): ${ratioChrNbLength}" \
    >> ${workdir}/5_coverage_anno/${SPECIES}.percent_Anno.txt
  echo "Number of Genes (after filtering): ${geneNumber}" >> ${workdir}/5_coverage_anno/${SPECIES}.percent_Anno.txt
  echo "Coverage of annotated Genes (bp, after filtering): ${geneCov}" >> ${workdir}/5_coverage_anno/${SPECIES}.percent_Anno.txt
  echo "Percentage of Genome Covered by Genes (%, after filtering): ${percentageGene}" >> ${workdir}/5_coverage_anno/${SPECIES}.percent_Anno.txt
  echo "Number of TEs: ${TEnumber}" >> ${workdir}/5_coverage_anno/${SPECIES}.percent_Anno.txt
  echo "Coverage of annotated TEs (bp): ${TECov}" >> ${workdir}/5_coverage_anno/${SPECIES}.percent_Anno.txt
  echo "Percentage of Genome Covered by TEs (%): ${percentageTE}" >> ${workdir}/5_coverage_anno/${SPECIES}.percent_Anno.txt

  # Output result in a tsv file
  echo -e "${SPECIES}\t${Gsize}\t${Gsize_wo500}\t${numberChr}\t${ratioChrNbLength}\t${geneNumber}\t${geneCov}\t${percentageGene}\t${TEnumber}\t${TECov}\t${percentageTE}" \
    >> ${OUTPUTFILE}

done < $INFILE



#########################################################################################################################

# 3. Create merging BED file for all species (for ridgeline plots)

TYPE_ANNO="gene TE"

for TYPE in $TYPE_ANNO; do 

  touch ${workdir}/allDarwin.merged_${TYPES}_cov.tsv

  # loop through each CSV file and append to merged_data.csv
  for FILE in ${workdir}/5_coverage_anno/*.${TYPE}_cov.bed; do
    
    SPECIES=$(basename "$FILE" | cut -d '.' -f 1)
    
    awk -v species="$SPECIES" 'NR==1{next} {print $0,"\t",species}' "$FILE" \
      >> ${workdir}/7_global_summaries/allDarwin.merged_${TYPE}_cov.tsv
  
  done

done


#############################################################################################################
# 4. Get the number of TE families per species

## 4.1. For families and subfamilies
while read -r _ SPECIES _ _; do
  
  awk '{print $8}' ${workdir}/2_BED_bedops/${SPECIES}.TE.sorted.bed | sort | uniq -c \
    > ${workdir}/6_TE_classes/${SPECIES}.count_TE.txt 

done < $INFILE	

## 4.2. For families (both count and cov_perc at the same time)

CATEGORIES="DNA	LINE LTR	PLE	RC	Retroposon	SINE	Uncertain	Unknown"

# Create summary file
touch ${workdir}/7_global_summaries/allDarwin.summary_TE_classes.txt
SUMMARY="${workdir}/7_global_summaries/allDarwin.summary_TE_classes.txt" 
> $SUMMARY # clear in case

# create header (5 columns)
echo -e "Species\tCategory\tCount\tCoverage\tCoverage_percent" >> ${SUMMARY}


while read -r _ SPECIES _ _; do

  # get genome size after trimming
  Gsize_wo500=$(awk '{sum += $2} END {print sum}' ${workdir}/0_inputs/fasta/${SPECIES}.wo500.fasta.fai)

  for category in $CATEGORIES ; do
    # calculate the number
    CATEGORY_COUNT=$(grep "${category}" ${workdir}/2_BED_bedops/${SPECIES}.TE.sorted.bed | wc -l) 

    if [ "$CATEGORY_COUNT" -eq "0" ]; then 
      
      echo "${category} is absent in ${SPECIES}"
      CATEGORY_COV="0" # fill the gaps
      CATEGORY_PERC="0" # fill the gaps
      
    else 
       
      # calculate the coverage (bp) 
      CATEGORY_COV=$(grep "${category}" ${workdir}/2_BED_bedops/${SPECIES}.TE.sorted.bed | awk '{sum += $3 - $2} END {print sum}')

      # calculate percentage
      CATEGORY_PERC=$(echo "scale=3; ($CATEGORY_COV * 100) / $Gsize_wo500" | bc)
    
    fi  
      
    # fill the summary by species
    echo -e "${SPECIES}\t${category}\t${CATEGORY_COUNT}\t${CATEGORY_COV}\t${CATEGORY_PERC}" \
      >> ${workdir}/6_TE_classes/${SPECIES}.TE_families.summary.txt

    # fill the global summary
    echo -e "${SPECIES}\t${category}\t${CATEGORY_COUNT}\t${CATEGORY_COV}\t${CATEGORY_PERC}" \
      >> ${SUMMARY}

  done

done < $INFILE