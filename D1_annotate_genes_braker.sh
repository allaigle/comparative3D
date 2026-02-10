#!/bin/bash
#
####################################################################################################################
#                                Annotate genes                    
#                   using BRAKER, bash and the SLURM scheduler                                                                             
#
# [GAOL] Annotate genome using BRAKER pipeline C (for reproducibility: rarely RNA-seq available for our species).
#
# [USAGE] sbatch -J D1 -c 8 --mem=40GB D1_annotate_genes_braker.sh SPECIES ASSEMBLY
#
# [Input] .wo500.fasta
# [Output] .gff3
#
# [Adaptation & History]
# Nov 2024 Alice Laigle (alice.laigle@gmail.com)
# May 2025 AL - added "PERL5LIB" and threads option.
# 
# NOTES:
#    1. Genomes are already downloaded (Cf. 'A1_doanload_genome_NCBI.sh') and
#       all sequences <500kb have been removed (Cf. 'A2_check_genome_quality_trim_wo500.sh')
#    2. BRAKER is already installed (their github: https://github.com/Gaius-Augustus/BRAKER)
#       and all needed tools and their PATH have been set up. 
#    3. fungal DB is already dowmloaded as proposed by BRAKER (see commented section 0)
###################################################################################################################
#
#SBATCH -c 8 # processor' number
#SBATCH --mem 40GB # memory

# Arguments
SPECIES=$1 # e.g, "Lflavum"
ASSEMBLY=$2 # e.g., "GCA_963580495.1"

source /data/alicel/miniconda3/etc/profile.d/conda.sh # to be changed
conda activate braker2

# Variables - TO BE CHANGED
basePATH="/data/alicel/chapter2"
workdir="${basePATH}/1_dataTreatment/3_annotations_braker"
export PATH=/data/alicel/miniconda3/bin:$PATH
export AUGUSTUS_CONFIG_PATH=/home/alicel/homebrew/Cellar/augustus/3.5.0_7/config
export BAMTOOLS_PATH=/home/alicel/homebrew/bin/bamtools 
export CDBTOOLS_PATH=/data/alicel/chapter2/4_tools/cdbfasta
export DIAMOND_PATH=/data/alicel/chapter2/4_tools
export GENEMARK_PATH=/data/alicel/chapter2/4_tools/GeneMark-ETP/bin/gmes
export PATH=/data/alicel/chapter2/4_tools/BRAKER/scripts:$PATH
export PATH=$PATH:/data/alicel/chapter2/4_tools/stringtie-2.2.3
export PATH=$PATH:/home/alicel/homebrew/bin/bedtools
export PATH=$PATH:/data/alicel/chapter2/4_tools/gffread
export PATH=$PATH:/data/alicel/chapter2/4_tools/compleasm_kit
export PROTHINT_PATH=/data/alicel/chapter2/4_tools/ProtHint/bin
export PYTHON3_PATH=/data/alicel/miniconda3/envs/braker2/bin/python3
export TSEBRA_PATH=/data/alicel/chapter2/4_tools/TSEBRA/bin
export MAKEHUB_PATH=/data/alicel/chapter2/4_tools/MakeHub
export PERL5LIB=/data/alicel/miniconda3/envs/braker2/lib/5.26.2

cd ${workdir}
mkdir -p ${SPECIES}

# 0. Fungal database
fungiDB="${workdir}/Fungi.fa"

# Add link of the FASTA file
ln -sf ${basePATH}/0_data/1_genomes/wo500/${SPECIES}_${ASSEMBLY}.wo500.fasta \
  ${workdir}/${SPECIES}_${ASSEMBLY}.wo500.fasta
FASTA="${SPECIES}_${ASSEMBLY}.wo500.fasta"

# Run BRAKER
braker.pl \
    --SPECIES=${SPECIES} \
    --genome=${FASTA} \
    --prot_seq=${fungiDB} \
    --gff3  \
    --workingdir=${workdir}/${SPECIES} \
    --threads=8 \
    --fungus

chmod -R 775 ${workdir}/${SPECIES}
