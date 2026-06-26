#!/usr/bin/env Rscript

# Goal: Plot Supplemental Figure  3 - Rigleline plot of gene/TE coverage
## Author: Alice Laigle

################################################################################

### Import libraries
library(ggplot2) # v4.0.0
library(ggridges) # v0.5.7
library(dplyr) # v1.1.4
library(patchwork) # v1.3.2
library(reshape2) # v1.4.4
theme_set(theme_bw()) # set nice theme
source("/Users/alicelaigle/PhD/fungalHiC/3_scripts/variables_names_colors_Fig1-3.R")

################################################################################

### Set PATHS
basePATH <- "/Users/alicelaigle/PhD/fungalHiC" # LOCAL

## path to phylo_tree_order_allDarwin.55.txt, df_info_perGenome.csv & categorized_models.txt
collectionPATH <- paste0(basePATH,"/0_data/0_collection/")
## path of summaries: TE classes and gene coverage
summAnnoPATH <- paste0(basePATH,"/2_analysis/1_annotations/7_global_summaries/")
## plot output path
plotPATH <- paste0(basePATH,"/5_final_results/SupplementalFigues/")

################################################################################

### Load files

df_phyloOrder <- read.table(paste0(collectionPATH, "phylo_speciesOrder_55.txt"),
                            sep = "\t", header= TRUE)
df_gene_mat <- read.table(paste0(summAnnoPATH,"allDarwin.merged_gene_cov.tsv"), 
                          sep = "\t")
df_TE_mat <- read.table(paste0(summAnnoPATH,"allDarwin.merged_TE_cov.tsv"), 
                        sep = "\t")

################################################################################

### 1. Ridgeline plot of gene coverage per 10kb bin 

# Rename % and species columns
names(df_gene_mat)[names(df_gene_mat) == "V7"] <- "percentage"
names(df_gene_mat)[names(df_gene_mat) == "V8"] <- "Species"

### Trim whitespaces first
df_gene_mat$Species <- trimws(df_gene_mat$Species)
df_phyloOrder$Species <- trimws(df_phyloOrder$Species)

# Remove the 3 species out of the threshold
df_gene_mat <- df_gene_mat %>% filter(!Species %in% dropped_species)

# Give full name of species
df_gene_mat <- df_gene_mat %>%
  left_join(df_phyloOrder %>% select(Species, Full_Species_Name), by = "Species")

# Reorder species
df_gene_mat$Full_Species_Name <- factor(df_gene_mat$Full_Species_Name, 
                                        levels = df_phyloOrder$Full_Species_Name)

### Put percentages as 0 to 100 instead of 0 to 1
df_gene_mat$percentage <- df_gene_mat$percentage * 100

### Plot

g_geneCov <- ggplot(df_gene_mat, 
                    aes(x = percentage, y = Full_Species_Name, fill = after_stat(x))) +
  geom_vline(xintercept = 20, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_density_ridges_gradient(alpha = 0.2, scale = 1, rel_min_height = 0.01) +
  scale_fill_viridis(name = "Percentage", option = "mako") +
  theme(axis.title = element_text(size = 12),
        axis.title.y = element_blank(), 
        axis.text.y = element_text(face = "italic", size = 10),
        axis.text.x = element_text(size = 9),
        strip.text.y.left = element_blank(),
        strip.text.x.top = element_text(color = "black"),
        strip.background.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        legend.position = "bottom") +
  xlab("Gene coverage/10kb bin (%)") +
  ylab("Species") 

g_geneCov

################################################################################

### 2. Ridgeline plot of TE coverage per 10kb bin 

# Rename % and species columns
names(df_TE_mat)[names(df_TE_mat) == "V7"] <- "percentage"
names(df_TE_mat)[names(df_TE_mat) == "V8"] <- "Species"

### Trim whitespaces first
df_TE_mat$Species <- trimws(df_TE_mat$Species)
df_phyloOrder$Species <- trimws(df_phyloOrder$Species)

# Remove the 3 species out of the threshold
df_TE_mat <- df_TE_mat %>% filter(!Species %in% dropped_species)

# Give full name of species
df_TE_mat <- df_TE_mat %>%
  left_join(df_phyloOrder %>% select(Species, Full_Species_Name), by = "Species")

# Reorder species
df_TE_mat$Full_Species_Name <- factor(df_TE_mat$Full_Species_Name, 
                                        levels = df_phyloOrder$Full_Species_Name)

### Put percentages as 0 to 100 instead of 0 to 1
df_TE_mat$percentage <- df_TE_mat$percentage * 100
### Plot

g_TEcov <- ggplot(df_TE_mat, 
                  aes(x = percentage, y = Full_Species_Name, fill = after_stat(x))) +
  geom_vline(xintercept = 20, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = 80, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_density_ridges_gradient(alpha = 0.2, scale = 1, rel_min_height = 0.01) +
  scale_fill_viridis(name = "Percentage", option = "mako") +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 9),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        strip.text.y.left = element_blank(),
        strip.text.x.top = element_text(color = "black"),
        strip.background.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        legend.position = "bottom") +
  xlab("TE coverage/10kb bin (%)")

g_TEcov

################################################################################

### Plot the two ridgelines plots together

nested <- (g_geneCov | g_TEcov) + 
  plot_layout(widths = c(1, 1)) 
nested

# save
ggsave(nested, 
       file=paste0(plotPATH,"SuppFigure3_ridgeline_coverage.png"), 
       width = 17, height = 25, units = "cm")
ggsave(nested, 
       file=paste0(plotPATH,"SuppFigure3_ridgeline_coverage.pdf"), 
       width = 17, height = 25, units = "cm")