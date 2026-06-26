#!/bin/R

# Goal: Plot the Supplemental Figure 5
## Phylogenetic ANOVAs between genetic features and 3D genome architectures

# Author: Alice Laigle
# Adapted from Camille Cornet's script
## Note: to be run before SuppFigure4_phyloreg.R - not sure if still relevant

library(ape) # for function corPagel
library(tidyverse) # v2.0.0
library(phytools) # v2.5.2
library(geiger) # v2.0.11
library(nlme) # for function gls - v3.1-168
library(ggplot2) # v4.0.0
library(ggpubr) # v0.6.1
library(patchwork) # v1.3.2
theme_set(theme_bw()) # set nice theme
source("/Users/alicelaigle/PhD/fungalHiC/3_scripts/variables_names_colors_Fig1-3.R")

################################################################################

### Set PATHS
basePATH <- "/Users/alicelaigle/PhD/fungalHiC/" # LOCAL

## path to df_info_genome.csv
collectionPATH <- paste0(basePATH,"0_data/0_collection/")
## path to ultrametric tree
treePATH <- paste0(basePATH,"1_dataTreatment/5_phyloTree/5_tree/")
## path of summary of types contacts 
contactsPATH <-  paste0(basePATH,"1_dataTreatment/5_conversions/")
## plot output path
plotPATH <- paste0(basePATH,"5_final_results/SupplementalFigues/")

################################################################################
### Load files

# tree
phylo_tree <- read.tree(paste0(treePATH,"Ultrametric_fungi_ordered"))
# genome information
df_info_genome <- read.delim(paste0(collectionPATH,"info_genome_TE_SuppFigures2_5.txt"),
                             sep = "\t", header = TRUE)
row.names(df_info_genome) <- df_info_genome$Species
# interactions - for QC
df_interactions <- read.table(paste0(contactsPATH, "allDarwin.summary_interactions.50000.corrected_ICE.txt"), 
                              sep="\t", header = TRUE)

################################################################################

# Check if same species in phylo and in dataset - for all plots
df_check <- name.check(phylo_tree, df_info_genome)
df_check
phylo_tree <- drop.tip(phylo_tree, df_check$tree_not_data)
#plotTree(phylo_tree)

###############################################################################

# Clean dataframes
## Add the Model column depending on the color code
df_info_genome <- df_info_genome %>%
  mutate(Model = ifelse(Model == "Bean", "Bean shape", Model))

df_info_genome <- df_info_genome %>%
  mutate(
    ColorCode = case_when(
      Model == "Rabl" ~ "#99CCCC",
      Model == "Bean shape" ~ "#08085E",
      Model == "Chromosome Territories" ~ "#990099",
      Model == "Rabl-Bean" ~ "#115B75", 
      Model == "Bean-Chromosome Territories" ~ "#560A7F",
      Model == "Chromosome Territories-Rabl" ~ "#9292b0",
    ))

# Reorder models
ordered_levels <- c("Rabl", 
                    "Bean shape", 
                    "Chromosome Territories",    
                    "Rabl-Bean", 
                    "Bean-Chromosome Territories", 
                    "Chromosome Territories-Rabl")

# Convert Model to a factor with specified levels
df_info_genome$Model <- factor(df_info_genome$Model, levels = ordered_levels)

## Add columns from genome info
df_info_genome <- merge(df_info_genome, 
                        df_interactions[, c("Species", "total_interactions")], 
                         by = "Species", 
                         all.x = TRUE) 



# Calculate cis/trans ratio
df_info_genome$Ratio <- (df_info_genome$Cis_Interactions_perc/df_info_genome$Trans_Interactions_perc)

# Make intermediates triangles
color_codes_to_reshape <- c("#115B75", "#560A7F", "#9292b0")
df_info_genome$ShapeType <- ifelse(df_info_genome$ColorCode %in% color_codes_to_reshape, 17, 16)

################################################################################
#### Phylogenetic ANOVAs

# 1. Is 3D model correlated to the ratio?

model3D <- setNames(df_info_genome$ColorCode, df_info_genome$Species)
ratio <- setNames(df_info_genome$Ratio, df_info_genome$Species)
model3D <- model3D[phylo_tree$tip.label]
ratio <- ratio[phylo_tree$tip.label]

t_model3D_ratio <- phylANOVA(phylo_tree, model3D, ratio, nsim=10000, posthoc=FALSE)
pval_ratio <- t_model3D_ratio$Pf # get Pr(>F) value # 0.822
pval_ratio <- round(pval_ratio, 3)
pval_ratio

plot_1 <- ggplot(data = df_info_genome, aes(x = Ratio, 
                                            y = Model, 
                                            color = Model)) +  
  geom_point(aes(size = 2.5, shape = as.factor(ShapeType))) +
  scale_color_manual(values = map_colorcode) + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 9),
        legend.position="none") + 
  labs(x = "Ratio of cis/trans interations", y = "3D model", tag = "A") + 
  annotate(geom="text", x=Inf, y = Inf,
           label=sprintf("italic('p')~'%s'",pval_ratio),
           parse=TRUE, hjust = 1.1, vjust = 1.5, size = 4, color = "black") +
  scale_y_discrete(limits=rev) 
plot_1

# No correlation between ratio and 3D models

################################################################################

## 2. Is the 3D model correlated to GC content?

GC <- setNames(df_info_genome$GC_percent_wo500_QUAST, df_info_genome$Species)
GC <- GC[phylo_tree$tip.label]
t_model3D_GC <- phylANOVA(phylo_tree, model3D, GC, nsim=10000, posthoc=FALSE)
pval_GC <- t_model3D_GC$Pf # get Pr(>F) value

plot_2 <- ggplot(data = df_info_genome, aes(x = GC_percent_wo500_QUAST,  
                                            y = Model, 
                                            color = Model)) +  
  geom_point(aes(size = 2.5, shape = as.factor(ShapeType))) +
  scale_color_manual(values = map_colorcode) + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 9),
        axis.title.y = element_blank(),
        axis.text.y= element_blank(),
        legend.position="none") + 
  labs(x = "GC content (%)", y = "3D model", tag = "B") + 
  annotate(geom="text", x=Inf, y = Inf,
           label=sprintf("italic('p')~'%s'",pval_GC),
           parse=TRUE, hjust = 1.1, vjust = 1.5, size = 4, color = "black") +
  scale_y_discrete(limits=rev) 
plot_2

# No correlation between GC content and 3D model 

################################################################################
## 3. Is the 3D model correlated to the chromosome number ?

n <- setNames(df_info_genome$Chromosome_Number, df_info_genome$Species)
n <- n[phylo_tree$tip.label]
t_model3D_nbChr <- phylANOVA(phylo_tree, model3D, n, nsim=10000, posthoc=FALSE)
pval_n <- t_model3D_nbChr$Pf # get Pr(>F) value

plot_3 <- ggplot(data = df_info_genome, aes(x = Chromosome_Number, 
                                            y = Model, 
                                            color = Model)) +  
  geom_point(aes(size = 2.5, shape = as.factor(ShapeType))) +
  scale_color_manual(values = map_colorcode) + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 9),
        axis.title.y = element_blank(),
        axis.text.y= element_blank(),
        legend.position="none") + 
  labs(x = "Chromosome number", y = "3D model", tag = "C") + 
  annotate(geom="text", x=Inf, y = Inf,
           label=sprintf("italic('p')~'%s'",pval_n),
           parse=TRUE, hjust = 1.1, vjust = 1.5, size = 4, color = "black") +
  scale_y_discrete(limits=rev) 
plot_3

# No correlation between chromosome number and 3D model 

################################################################################

## 4. Is the 3D model correlated to gene coverage ?

gene <- setNames(df_info_genome$Gene_cov_perc, df_info_genome$Species)
gene <- gene[phylo_tree$tip.label]
t_model3D_gene <- phylANOVA(phylo_tree, model3D, gene, nsim=10000, posthoc=FALSE)
pval_gene <- t_model3D_gene$Pf # get Pr(>F) value

plot_4 <- ggplot(data = df_info_genome, aes(x = Gene_cov_perc, 
                                            y = Model, 
                                            color = Model)) +  
  geom_point(aes(size = 2.5, shape = as.factor(ShapeType))) +
  scale_color_manual(values = map_colorcode) + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 9),
        legend.position="none") + 
  labs(x = "Gene content (%)", y = "3D model", tag = "D") + 
  annotate(geom="text", x=Inf, y = Inf,
           label=sprintf("italic('p')~'%s'",pval_gene),
           parse=TRUE, hjust = 1.1, vjust = 1.5, size = 4, color = "black") +
  scale_y_discrete(limits=rev) 
plot_4

# No correlation between gene coverage and 3D model 

################################################################################

## 5. Is the genome size correlated to specific TEs?

### LTR ########################################################################

LTR <- setNames(df_info_genome$LTR, df_info_genome$Species)
LTR <- LTR[phylo_tree$tip.label]
t_model3D_LTR <- phylANOVA(phylo_tree, model3D, LTR, nsim=10000, posthoc=FALSE)
pval_LTR <- t_model3D_LTR$Pf # get Pr(>F) value

plot_5 <- ggplot(data = df_info_genome, aes(x = LTR, 
                                            y = Model, 
                                            color = Model)) +  
  geom_point(aes(size = 2.5, shape = as.factor(ShapeType))) +
  scale_color_manual(values = map_colorcode) + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 9),
        axis.title.y = element_blank(),
        axis.text.y= element_blank(),
        legend.position="none") + 
  labs(x = "LTR content (%)", y = "3D model", tag = "E") + 
  annotate(geom="text", x=Inf, y = Inf,
           label=sprintf("italic('p')~'%s'",pval_LTR),
           parse=TRUE, hjust = 1.1, vjust = 1.5, size = 4, color = "black") +
  scale_y_discrete(limits=rev) 
plot_5

# No correlation between LTR and 3D models

### DNA ########################################################################

DNA <- setNames(df_info_genome$DNA, df_info_genome$Species)
DNA <- DNA[phylo_tree$tip.label]
t_model3D_DNA <- phylANOVA(phylo_tree, model3D, DNA, nsim=10000, posthoc=FALSE)
pval_DNA <- t_model3D_DNA$Pf # get Pr(>F) value

plot_6 <- ggplot(data = df_info_genome, aes(x = DNA, 
                                            y = Model, 
                                            color = Model)) +  
  geom_point(aes(size = 2.5, shape = as.factor(ShapeType))) +
  scale_color_manual(values = map_colorcode) + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 9),
        axis.title.y = element_blank(),
        axis.text.y= element_blank(),
        legend.position="none") + 
  labs(x = "DNA content (%)", y = "3D model", tag = "F") + 
  annotate(geom="text", x=Inf, y = Inf,
           label=sprintf("italic('p')~'%s'",pval_DNA),
           parse=TRUE, hjust = 1.1, vjust = 1.5, size = 4, color = "black") +
  scale_y_discrete(limits=rev) 
plot_6

# Correlation between DNA and 3D models

### LINE ########################################################################

LINE <- setNames(df_info_genome$LINE, df_info_genome$Species)
LINE <- LINE[phylo_tree$tip.label]
t_model3D_LINE <- phylANOVA(phylo_tree, model3D, LINE, nsim=10000, posthoc=FALSE)
pval_LINE <- t_model3D_LINE$Pf # get Pr(>F) value

plot_7 <- ggplot(data = df_info_genome, aes(x = LINE, 
                                            y = Model, 
                                            color = Model)) +  
  geom_point(aes(size = 2.5, shape = as.factor(ShapeType))) +
  scale_color_manual(values = map_colorcode) + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 9),
        legend.position="none") + 
  labs(x = "LINE content (%)", y = "3D model", tag = "G") + 
  annotate(geom="text", x=Inf, y = Inf,
           label=sprintf("italic('p')~'%s'",pval_LINE),
           parse=TRUE, hjust = 1.1, vjust = 1.5, size = 4, color = "black") +
  scale_y_discrete(limits=rev) 
plot_7

# No correlation between LINE and 3D models

### PLE ########################################################################

PLE <- setNames(df_info_genome$PLE, df_info_genome$Species)
PLE <- PLE[phylo_tree$tip.label]
t_model3D_PLE <- phylANOVA(phylo_tree, model3D, PLE, nsim=10000, posthoc=FALSE)
pval_PLE <- t_model3D_PLE$Pf # get Pr(>F) value

plot_8 <- ggplot(data = df_info_genome, aes(x = PLE, 
                                            y = Model, 
                                            color = Model)) +  
  geom_point(aes(size = 2.5, shape = as.factor(ShapeType))) +
  scale_color_manual(values = map_colorcode) + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 9),
        axis.title.y = element_blank(),
        axis.text.y= element_blank(),
        legend.position="none") + 
  labs(x = "PLE content (%)", y = "3D model", tag = "H") + 
  annotate(geom="text", x=Inf, y = Inf,
           label=sprintf("italic('p')~'%s'",pval_PLE),
           parse=TRUE, hjust = 1.1, vjust = 1.5, size = 4, color = "black") +
  scale_y_discrete(limits=rev) 
plot_8

# No correlation between PLE and 3D models

### RC ########################################################################

RC <- setNames(df_info_genome$RC, df_info_genome$Species)
RC <- RC[phylo_tree$tip.label]
t_model3D_RC <- phylANOVA(phylo_tree, model3D, RC, nsim=10000, posthoc=FALSE)
pval_RC <- t_model3D_RC$Pf # get Pr(>F) value

plot_9 <- ggplot(data = df_info_genome, aes(x = RC, 
                                            y = Model, 
                                            color = Model)) +  
  geom_point(aes(size = 2.5, shape = as.factor(ShapeType))) +
  scale_color_manual(values = map_colorcode) + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 9),
        axis.title.y = element_blank(),
        axis.text.y= element_blank(),
        legend.position="none") + 
  labs(x = "RC content (%)", y = "3D model", tag = "I") + 
  annotate(geom="text", x=Inf, y = Inf,
           label=sprintf("italic('p')~'%s'",pval_RC),
           parse=TRUE, hjust = 1.1, vjust = 1.5, size = 4, color = "black") +
  scale_y_discrete(limits=rev) 
plot_9

# No correlation between RC and 3D models

### SINE ########################################################################

SINE <- setNames(df_info_genome$SINE, df_info_genome$Species)
SINE <- SINE[phylo_tree$tip.label]
t_model3D_SINE <- phylANOVA(phylo_tree, model3D, SINE, nsim=10000, posthoc=FALSE)
pval_SINE <- t_model3D_SINE$Pf # get Pr(>F) value

plot_10 <- ggplot(data = df_info_genome, aes(x = SINE, 
                                            y = Model, 
                                            color = Model)) +  
  geom_point(aes(size = 2.5, shape = as.factor(ShapeType))) +
  scale_color_manual(values = map_colorcode) + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 9),
        legend.position="none") + 
  labs(x = "SINE content (%)", y = "3D model", tag = "J") + 
  annotate(geom="text", x=Inf, y = Inf,
           label=sprintf("italic('p')~'%s'",pval_SINE),
           parse=TRUE, hjust = 1.1, vjust = 1.5, size = 4, color = "black") +
  scale_y_discrete(limits=rev) 
plot_10

# No correlation between SINE and 3D models

################################################################################
plot_blank <- ggplot() + theme_void() # create blanck plot
nested <- (plot_1 | plot_2 | plot_3) /
  (plot_4 | plot_5 | plot_6) /
  (plot_7 | plot_8 | plot_9) /
  (plot_10 | plot_blank | plot_blank) +
  plot_layout(guides = "collect") & theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
nested

# Save
ggsave(nested, 
       file=paste0(plotPATH,"SuppFigure5_phyloANOVA.pdf"), 
       width = 20, height = 20, units = "cm")

plot_1 <- plot_1 + theme(axis.text.y = element_blank())
plot_4 <- plot_4 + theme(axis.text.y = element_blank())
plot_7 <- plot_7 + theme(axis.text.y = element_blank())
plot_10 <- plot_10 + theme(axis.text.y = element_blank())
nested <- (plot_1 | plot_2 | plot_3) /
  (plot_4 | plot_5 | plot_6) /
  (plot_7 | plot_8 | plot_9) /
  (plot_10 | plot_blank | plot_blank) +
  plot_layout(guides = "collect") & theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
ggsave(nested, 
       file=paste0(plotPATH,"SuppFigure5_phyloANOVA.wo_ytext.pdf"), 
       width = 15, height = 22.5, units = "cm")
