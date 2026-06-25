#!/usr/bin/env Rscript

# Goal: Plot panel 1 - Plot phylogenetic tree
## Author: Alice Laigle

################################################################################
# Import libraries
library(readr) # v2.1.5
library(ape) # v5.8.1
library(phytools) # v2.5.2
library(ggplot2) # v4.0.0
library(ggtree) # v3.16.3 
library(dplyr) # v1.1.4
source("/Users/alicelaigle/PhD/fungalHiC/3_scripts/variables_names_colors_Fig1-3.R")

################################################################################

### Set paths

# LOCAL
## path of tree file
inTreePATH <- "/Users/alicelaigle/PhD/fungalHiC/1_dataTreatment/5_phyloTree/5_tree/"
## plot output path
plotPATH <- "/Users/alicelaigle/PhD/fungalHiC/5_final_results/Figure1_tree_models/"

################################################################################

### Load tree
chr_tree <- read.tree(paste0(inTreePATH,"Ultrametric_fungi_ordered"))

################################################################################
### Prepare data to plot tree

# Replace the tip labels in the tree object
chr_tree$tip.label <- name_mapping[chr_tree$tip.label] # name_mapping is sourced

# Remove 5 species (those under 1 million valid pairs wihtout duplicates)
chr_tree_cleaned <- drop.tip(chr_tree, dropped_tip_FullName) # dropped_tip_FullName is sourced

################################################################################

### Define palettes and phyla

### Palette
palettePhylum <- c("Zoopagomycota" = "black", "Mucoromycota" = "orange", 
                   "Ascomycota" = "#6BBDFF", "Basidiomycota" = "#cc4778")
# Add tip colors 
tip_colors <- palettePhylum[phyla[chr_tree_cleaned$tip.label]] # phyla is sourced

# Make it as simmap
set.seed(12) 
simmap_tree <- make.simmap(chr_tree_cleaned, phyla, nsim=1)

################################################################################
### Plot tree as linear
pdf(file = paste0(plotPATH,"Figure1A_tree_linear_58species.pdf"), 
    width = 5, height = 10)

plotSimmap(simmap_tree, ftype="i", 
           colors=tip_colors, offset = 0.5) +
  tiplabels(pch = map_shapes, col = NA, bg = map_colors, cex = 1.75)

dev.off()