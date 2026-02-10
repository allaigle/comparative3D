This repository contains all the scripts used in the study "Comparative 3D genome architecture analyses reveal rapid turnover and the emergence of transitory stages in the fungal kingdom". It uses the SLURM scheduler and contains scripts needing variable changes.

# Structure of the project 

For clarification, gave the most important structural folders.

```{bash}
$basePATH 
├── 0_data
│   ├── 0_collection
│   ├── 1_genomes
├── 1_dataTreatment
│   ├── 0_QUAST
│   ├── 1_3DGB
│   ├── 2_phyloTree
│   ├── 3_annotations_braker
│   ├── 4_annotations_EarlGrey
│   ├── 5_conversions
│   ├── 6_stats_HiC_Pro
├── 2_analysis
│   ├── 1_annotations
│   ├── 2_cross_hic_anno
├── 3_scripts
└── 4_tools
```

# A - Collection of genomes

Public fungal genomes are mostly from [The Darwin Tree of Life (DToL) Project Consortium, 2022](https://www.pnas.org/doi/full/10.1073/pnas.2115642118). The GenBank accessions used in this study can be found in the Supplemental Table 1.

## Download genomes

Requires: `datasets` from [NCBI](https://github.com/ncbi/datasets).

Input: `allDarwin_phylum_species_assembly_SRA.filtered_60.txt`, TXT file containing at least 3 columns with PHYLUM (abbreviated, *e.g.*, basidio), SPECIES (abbreviated, *e.g.*, Abisporus) and ASSEMBLY (GenBank ID).

```{bash}
bash A1_download_datasets_NCBI.sh # as an array
```

Note: for the three outgroup species, used directly the sbatch command inside it as

```{bash}
# Example:
basePATH="/data/alicel/chapter2"
datasets="${basePATH}/4_tools/datasets"
PHYLUM="zoopago" # abbreviation of Zoopagomycota
SPECIES="Cmojavensis" # abbreviation of Coemansia mojavensis
ASSEMBLY="GCA_025766245.1" # GenBank assembly

cd ${basePATH}/0_data/1_genomes/slurm

sbatch -J A1_${SPECIES} -c 1 --mem=8MB --wrap="mkdir -p ../${PHYLUM}/${SPECIES} ;\
     $datasets download genome accession $ASSEMBLY --include genome,seq-report \
     --filename ../${PHYLUM}/${SPECIES}/${SPECIES}_${ASSEMBLY}_dataset.zip ;\
     cd ../${PHYLUM}/${SPECIES} ;\
     unzip ${SPECIES}_${ASSEMBLY}_dataset.zip"
```

## Genome quality

Removing sequences smaller than 500kb (for 3DGB later, not working if too small contigs/scaffolds) requires SAMtools [Danecek *et al.,* 2021](https://academic.oup.com/gigascience/article/10/2/giab008/6137722?login=false) and [seqtk](https://github.com/lh3/seqtk).

Input: `allDarwin_phylum_species_assembly_SRA.filtered_60.txt`.

```{bash}
sbatch -J A2 -c 1 --mem=16MB A2_check_genome_quality_trim_wo500.sh
```


Quality control using QUAST ([Gurevich *et al.,* 2013](https://academic.oup.com/bioinformatics/article/29/8/1072/228832)) and miniconda (used alias to conda).

Input: `allDarwin_phylum_species_assembly_SRA.filtered_60.txt`.

```{bash}
sbatch -J A3 -c 1 --mem=16MB A3_QUAST_allDarwin.sh
```

Note: For this part, did the same with `zoopagomycota_species_assembly.txt` as input to get the outgroups for the generation of the phylogenetic tree (*Cf.* Construct phylogenetic tree).


# B - Collection and process of Hi-C data using 3DGB pipeline

<img align="right" width="200px" 
    src="Bsubappendiculatus_10kb.gif"
    alt="3D structure of the Butyriboletus_subappendiculatus genome at 10 kb resolution">

3D Genome Builder (3DGB) is a snakemake workflow developped by [Poisignon *et al.* (2022)](https://doi.org/10.1093/nargab/lqad104), treating Hi-C data (requires only FASTA and FASTQ.GZ files), controlling their quality and generating PDB file formats. The PDB can then be uploaded into different softwares to visualize the 3D models. In our case, we used [Mol* Viewer](https://molstar.org/viewer/) (Element Index option in Spacefill category for the colors, giving a gradient where the longest chromosome is red and the smallest blue).


**Important note:** the tree generation and the annotations have been done in the same time since 3DGB took some time to be set up and run on our cluster, espacially with this number of species.

## Prepare workdir 

Requires: `prefetch` and `fasterq-dump` from [SRA Toolkit](https://github.com/ncbi/sra-tools).

```{bash}
# [USAGE] sbatch -J B1 -c 1 --mem=2MB B1_prepare_3DGB.sh PHYLUM SPECIES \
## FULL_SPECIES_NAME ASSEMBLY SRA REFULL_SPECIES_NAME ASSEMBLY SRA RE
# Example: 
sbatch -J B1 -c 1 --mem=2MB B1_prepare_3DGB.sh muco Eparvispora \
  "Entomortierella parvispora" GCA_963556305.1 ERR11577535 "^GATC G^ANTC C^TNAG T^TAA" or "^GATC" 
```

Notes: 
1) the 3DGB configuration file is "pre"-modified for an automation purpose (*Cf.* `config_template_fungalHiC.yml`),
2) `config_template_fungalHiC.yml` is actually in a `configYML` folder, and
3) if we have digged into the different resolutions, at the end, we only used the 10 and 50kb resolutions.

## Run and clean workdir 

```{bash}
#[USAGE] bash B2_run_3DGB.sh PHYLUM SPECIES SRA
# Example: 
bash B2_run_3DGB.sh muco Eparvispora ERR11577535 
```

Note: `B2_run_3DGB.sh` contains an sbatch dependency, including `B3_clean_3DGB_workdir.sh`, which is very specific to our cluster since we had to remove heavy files as soon as possible.


# C - Generation of the phylogenetic tree

Scripts of this section are based on Toby Baril's ones. 

After finding three fungal genomes of good quality out of the three phyla present in the DToL, in our case from the Zoopagomycota phylum, we could run BUSCO for the 63 species.

## Step 1: Prepare and run BUSCO

Requires BUSCO (v5.8.0; [Manni *et al.,* 2021](https://academic.oup.com/mbe/article/38/10/4647/6329644)).

Input: `phylotree_phylum_species_assembly.63.txt`.

```{bash}
sbatch -J C1 -c 1 --mem=2MB C1_prepare_run_BUSCO.sh
```

Note: `C1_prepare_run_BUSCO.sh` contains sbatch dependencies and `C2_run_busco_TBaril.sh`.

## Step 2: Extract and filter complete BUSCO genes from runs

Briefly: 
1) Extract any gene from each directory that is complete,
2) sort the complete IDs and count how many genomes they appeared in, and
3) keep only those present in at least 3 genomes.

Input: `phylotree_phylum_species_assembly.63.txt`.

```{bash}
sbatch -J C3 -c 1 --mem=16MB C3_extract_filter_BUSCO_TBaril.sh
```

## Step 3: Prepare files for alignements

Briefly: 
1) Copy amino acid sequences and name files with species name in them,
2) add species names to FASTA headers, and
3) generate a single file for each BUSCO ID.

Input: `phylotree_phylum_species_assembly.63.txt`.

```{bash}
sbatch -J C4 -c 1 --mem=16MB C4_prepare_alignement_TBaril.sh
```

## Step 3.5: Check point and deletion of individual FASTA files

- Check the fasta headers are consistent (this is really important for later).
- Check the general file structure looks okay.
-> **Done**.

Note: individual FASTA files are deleted as these are still saved in the BUSCO directories.

```{bash}
cd /data/alicel/chapter2/1_dataTreatment/2_phyloTree/2_busco_aa

for file in *.faa ;
	do rm ${file} ; 
done
```

## Step 4: Align peptide sequences for each gene

Requires Mafft (v7.626; [Katoh and Standley, 2013](https://academic.oup.com/mbe/article/30/4/772/1073398)), in an environement.

```{bash}
sbatch -J C5 -c 1 --mem=4MB C5_prepare_run_mafft_TBaril.sh
```

Note: `C5_prepare_run_mafft_TBaril.sh` contains `C6_run_mafft_TBaril.sh`.

Output: `${basePATH}/1_dataTreatment/2_phyloTree/complete_aligned_ids.txt`.

## Step 5: Edit headers

The goal here is to have headers being homogeneous in every alignment file.

```{bash}
sbatch -J C7 -c 1 --mem=4MB C7_edit_headers_TBaril.sh
```

## Step 6: Create a supermatrix

Requires PhyKIT (v2.0.1; [Steenwyk *et al.,* 2021](https://academic.oup.com/bioinformatics/article/37/16/2325/6131675))

Briefly: 
1) create a list of alignment files,
2) concatenate the alignments, and
3) trim uninformative regions.

```{bash}
sbatch -J C8 -c 1 --mem=32MB C8_create_supermatrix_TBaril.sh
```

## Step 7: Run the tree with bootstraps

Requires IQ-Tree2 ([Minh *et al.,* 2020](https://academic.oup.com/mbe/article/37/5/1530/5721363?login=false)).

Input: `fungi_supermatrix.fa.clipkit`.

```{bash}
sbatch -J C9 -c 1 --mem=2MB C9_run_iqtree2_TBaril.sh
```

Notes:
- contains a sbatch command for running IQ-Tree2 with much more ressources, 
- takes the three Zoopagomycota species as outgroups.


# D - Gene and TE annotations

## Gene annotation

We used the pipeline C of BRAKER2 ([Brůna *et al.,* 2021](https://doi.org/10.1093/nargab/lqaa108)) to all species for a reproducibility purpose since only few species were having RNA-seq data available. The fungal database we used 

```{bash}
# [USAGE] sbatch -J D1 -c 8 --mem=40GB D1_annotate_genes_braker.sh SPECIES ASSEMBLY
# Example
sbatch -J D1 -c 8 --mem=40GB D1_annotate_genes_braker.sh Lflavum GCA_963580495.1
```

## TE annotation
