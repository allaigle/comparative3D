#!/bin/bash
#
####################################################################################################################
#                          Prepare folders, files and config for 3DGB                        
#                        using SRA Toolkit, bash and the SLURM scheduler                                                                             
#
# [USAGE] bash B1_prepare_3DGB.sh PHYLUM SPECIES FULL_SPECIES_NAME ASSEMBLY SRA RE
#
# [Adaptation & History]
# Sept-Oct 2024 - Alice Laigle (alice.laigle@unine.ch)
# Nov 2024 - AL - adapted to genomes without <500kb & created dependencies
# May 2025 - AL - Adapt to the new cluster system
#
# [Input] ASSEMBLY and SRA accessions
# [Output] FASTA, FASTQ.QZ and YML config file
# 
# NOTES:
#    1. To have as much reproducibility as possible, use FASTA without <500kb chr/contigs/scaffolds for all. 
#    2. 3DGB, datasets & SRA Toolkit are already installed
#    3. Genome is already downloaded.
#    4. 3DGB config file is already modified, see 'config_template_fungalHiC.yml'. 
###################################################################################################################
#
#SBATCH -c 1 # processor' number
#SBATCH --mem 2MB # memory

# Arguments
PHYLUM=$1 # e.g., "muco"
SPECIES=$2 # e.g., "Eparvispora"
FULL_SPECIES_NAME=$3 # e.g., "Entomortierella parvispora"
ASSEMBLY=$4 # e.g.,  "GCA_963556305.1"
SRA=$5 # e.g., "ERR11577535"
RE=$6 # e.g., "^GATC G^ANTC C^TNAG T^TAA" or "^GATC"

# Variable - TO BE CHANGED
basePATH="/data/alicel/chapter2"
workdir="${basePATH}/1_dataTreatment/1_3DGB"
fasterq-dump="/home/alicel/tools/fasterq-dump"

# 1. Download FASTQ files in 3DGB workdir folder
mkdir -p ${workdir}/${SPECIES}/fastq_files 
cd ${workdir}/${SPECIES}/fastq_files

sbatch -J prefetch_${SPECIES} -c 1 --mem 16MB \
  --wrap="cd ${workdir}/${SPECIES}/fastq_files ;\
  prefetch ${SRA} --max-size 10000000000"

# dependency 1
sbatch -J fasterq-dump_${SPECIES} -c 1 --mem 2GB \
  --dependency=$(squeue --noheader --format %i --name prefetch_${SPECIES}) \
  --wrap="cd ${workdir}/${SPECIES}/fastq_files/${SRA} ;\
   ${fasterq-dump} --split-files --include-technical ${SRA}"

# 2. Format FASTQ files as required by 3DGB
# dependency 2
sbatch -J mv_clean_${SPECIES} -c 1 --mem 8MB \
  --dependency=$(squeue --noheader --format %i --name fasterq-dump_${SPECIES}) \
  --wrap="cd ${workdir}/${SPECIES}/fastq_files/${SRA} ;\
  mv ${SRA}_1.fastq ${SRA}_R1.fastq ;\
  mv ${SRA}_2.fastq ${SRA}_R2.fastq ;\
  rm ${SRA}.sra"

# dependency 3
sbatch -J gzipFASTQ_1_${SPECIES} -c 1 --mem 8MB \
  --dependency=$(squeue --noheader --format %i --name mv_clean_${SPECIES}) \
  --wrap="cd ${workdir}/${SPECIES}/fastq_files/${SRA} ;\
   gzip ${SRA}_R1.fastq"

sbatch -J gzipFASTQ_2_${SPECIES} -c 1 --mem 8MB \
  --dependency=$(squeue --noheader --format %i --name mv_clean_${SPECIES}) \
  --wrap="cd ${workdir}/${SPECIES}/fastq_files/${SRA} ;\
   gzip ${SRA}_R2.fastq"

# 3. Get genome in 3DGB workdir folder
cd ${workdir}/${SPECIES}
cp ${basePATH}/0_data/1_genomes/wo500/${SPECIES}_${ASSEMBLY}.wo500.fasta genome.fasta
chmod 775 genome.fasta 

# 4. Adapt 3DGB config file (here for Darwin subset, no need of changing RE)
cp ${workdir}/configYML/config_template_fungalHiC.yml \
  ${workdir}/configYML/${PHYLUM}/${SPECIES}_config.yml
sed -i 's/workdir/'${SPECIES}'/g' ${workdir}/configYML/${PHYLUM}/${SPECIES}_config.yml
sed -i 's/FULLNAME/'"${FULL_SPECIES_NAME}"'/g' ${workdir}/configYML/${PHYLUM}/${SPECIES}_config.yml
sed -i 's/ACCESSION/'${SRA}'/g' ${workdir}/configYML/${PHYLUM}/${SPECIES}_config.yml
sed -i 's/RESTRICTION_ENZYMES/'"${RE}"'/g' ${workdir}/configYML/${PHYLUM}/${SPECIES}_config.yml
cat ${workdir}/configYML/${PHYLUM}/${SPECIES}_config.yml 

chmod -R 775 ${workdir}/${SPECIES}/*