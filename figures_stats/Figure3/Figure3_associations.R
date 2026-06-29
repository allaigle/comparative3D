#!/usr/bin/env Rscript

# Goal: Plot the Figure 3 and get stats - phylogenetic regressions

# Author: Alice Laigle
# Adapted from Camille Cornet's script

library(ape) # for function corPagel - v5.8-1
library(tidyverse) # v2.0.0
library(phytools) # v2.5.2
library(geiger) # v2.0.11
library(nlme) # for function gls - v3.1-168
library(ggplot2) # v4.0.0
library(ggpubr) # v0.6.1
theme_set(theme_bw()) # set nice theme
source("/Users/alicelaigle/PhD/fungalHiC/3_scripts/variables_names_colors_Fig1-3.R")

################################################################################

### Set PATHS
basePATH <- "/Users/alicelaigle/PhD/fungalHiC/" # LOCAL

## path to phylo_tree_order_allDarwin.55.txt, df_info_genome.csv & categorized_models.txt
collectionPATH <- paste0(basePATH,"0_data/0_collection/")
## path to ultrametric tree
treePATH <- paste0(basePATH,"1_dataTreatment/5_phyloTree/5_tree/")
## plot output path
plotPATH <- paste0(basePATH,"5_final_results/Figure3_associations/")

################################################################################

### Load files

# tree
phylo_tree <- read.tree(paste0(treePATH,"Ultrametric_fungi_ordered"))

# genome information
df_info_genome <- read.delim(paste0(collectionPATH,"info_genome.55species.Figure3.csv"),
                             sep = ",", header = TRUE, row.names = 1)
df_info_genome$Species <- row.names(df_info_genome)

################################################################################

# Check if same species in phylo and in dataset
df_check <- name.check(phylo_tree, df_info_genome)
df_check
phylo_tree <- drop.tip(phylo_tree, df_check$tree_not_data)

# Define a variance-covariance structure
pagel <- corPagel(1, phy = phylo_tree, fixed = FALSE, form = ~Species)

################################################################################

### Prepare data
# Add the ColorCode column depending on the Model to df_info_genome
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

# Calculate cis/trans ratio
df_info_genome$Ratio <- (df_info_genome$Cis_Interactions_perc/df_info_genome$Trans_Interactions_perc)

# Make intermediates triangles
color_codes_to_reshape <- c("#115B75", "#560A7F", "#9292b0")
df_info_genome$ShapeType <- ifelse(df_info_genome$ColorCode %in% color_codes_to_reshape, 17, 16)


################################################################################
#### Phylogenetic Generalized Least Squares

## 1. Is the genome size correlated to TE content ?

# run model
pagel <- corPagel(1, phy = phylo_tree, fixed = FALSE, form = ~Species)
model_gsize_TE <- gls(GenomeSize_bp_wo500 ~ TE_cov_perc,
                      data = df_info_genome, association = pagel)
summary(model_gsize_TE)

# extract values
p_value_size_TE <- summary(model_gsize_TE)$tTable[2, "p-value"] # get p-value
p_value_size_TE <- round(p_value_size_TE, 11) # write in shorter way for plotting

aic_value_size_TE <- AIC(model_gsize_TE) # get AIC
aic_value_size_TE
bic_value_size_TE <- BIC(model_gsize_TE) # get  BIC
bic_value_size_TE

# Get genome size into Mbp instead of bp
df_info_genome$GenomeSize_Mbp <- (df_info_genome$GenomeSize_bp_wo500)/1000000

# Reorder models
ordered_levels <- c("Rabl", 
                    "Bean shape", 
                    "Chromosome Territories",    
                    "Rabl-Bean", 
                    "Bean-Chromosome Territories", 
                    "Chromosome Territories-Rabl")

# Convert Model to a factor with specified levels
df_info_genome$Model <- factor(df_info_genome$Model, levels = ordered_levels)

# plot
plot_1 <- ggplot(data = df_info_genome, aes(x = GenomeSize_Mbp, 
                                            y = TE_cov_perc,
                                            colour = Model)) +  
  geom_point(aes(size = 2.5, shape = as.factor(ShapeType))) +
  scale_color_manual(values = map_colorcode) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 9),
        legend.position = "none"
        ) + 
  labs(x = "Genome size (Mbp)", y = "TE content (%)", color = "3D model",
       tag = "A") + 
  annotate(geom="text", x = 150, y = 15,
           label=sprintf("italic('p')~'%s'",p_value_size_TE),
           parse=TRUE, hjust = 1.1, vjust = 1.5, size = 4, color = "black") +
  annotate("text", x = 160, y = 10, 
           label = paste0("AIC = ",round(aic_value_size_TE, 2)), 
           hjust = 1.1, vjust = 1.5, size = 4, color = "black") +
  annotate("text", x = 160, y = 5, 
           label = paste0("BIC = ",round(bic_value_size_TE, 2)), 
           hjust = 1.1, vjust = 1.5, size = 4, color = "black")
  
plot_1

# Very strong association between the genome size and the TE content

################################################################################

## 2. Is the ratio of cis/trans correlated with TE content ?

model_ratio_TEcov <- gls(Ratio ~ TE_cov_perc,
                         data = df_info_genome, association = pagel)
summary(model_ratio_TEcov)

# extract values
p_value_ratio_TEcov <- summary(model_ratio_TEcov)$tTable[2, "p-value"]
p_value_ratio_TEcov <- round(p_value_ratio_TEcov, 4) # checked manually
aic_value_ratio_TEcov <- AIC(model_ratio_TEcov) # get AIC
bic_value_ratio_TEcov <- BIC(model_ratio_TEcov) # get  BIC


plot_2 <- ggplot(data = df_info_genome, aes(x = Ratio, 
                                              y = TE_cov_perc, 
                                            colour = Model)) +  
  geom_point(aes(size = 2.5, shape = as.factor(ShapeType))) +
  scale_color_manual(values = map_colorcode) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 9)) + 
  labs(x = "Ratio cis/trans interactions", y = "TE content (%)", 
       tag = "B", color = "3D model") + 
  annotate(geom="text", x=11,y=15,
           label=sprintf("italic('p')~'%s'",p_value_ratio_TEcov),
           parse=TRUE, hjust = 1.1, vjust = 1.5, size = 4, color = "black") +
  annotate("text", x = 12, y = 10, 
           label = paste0("AIC = ",round(aic_value_ratio_TEcov, 2)), 
           hjust = 1.1, vjust = 1.5, size = 4, color = "black") +
  annotate("text", x = 12, y = 5, 
           label = paste0("BIC = ",round(bic_value_ratio_TEcov, 2)), 
           hjust = 1.1, vjust = 1.5, size = 4, color = "black")
plot_2
# No association between the ratio cis/trans and TE content

################################################################################
#### Phylogenetic ANOVAs

## 3. Is 3D model correlated to genome size?
### note: ColorCode refers to a 3D model as in Figure 1

model3D <- setNames(df_info_genome$ColorCode, df_info_genome$Species)
size <- setNames(df_info_genome$GenomeSize_bp_wo500, df_info_genome$Species)
model3D <- model3D[phylo_tree$tip.label]
size <- size[phylo_tree$tip.label]

t_model3D_size <- phylANOVA(phylo_tree, model3D, size, nsim=10000, posthoc=FALSE)
p_value_size <- t_model3D_size$Pf # get Pr(>F) value
p_value_size <- round(p_value_size, 3) # checked manually

plot_3 <- ggplot(data = df_info_genome, aes(x = GenomeSize_Mbp, 
                                            y = Model, 
                                            color = Model)) +  
  geom_point(aes(size = 2.5, shape = as.factor(ShapeType))) +
  scale_color_manual(values = map_colorcode) + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 9),
        legend.position="none") + 
  labs(x = "Genome size (Mbp)", y = "3D model", tag = "C") + 
  annotate(geom="text", x=Inf, y = Inf,
           label=sprintf("italic('p')~'%s'",p_value_size),
           parse=TRUE, hjust = 1.1, vjust = 1.5, size = 4, color = "black") +
  scale_y_discrete(limits=rev) 
plot_3

# No association between genome size and 3D models

################################################################################

## 4. Is the 3D model correlated to TE content ?

TE <- setNames(df_info_genome$TE_cov_perc, df_info_genome$Species)
TE <- TE[phylo_tree$tip.label]
t_model3D_TE <- phylANOVA(phylo_tree, model3D, TE, nsim=10000, posthoc=FALSE)
p_value_TE <- t_model3D_TE$Pf # get Pr(>F) value
p_value_TE <- round(p_value_TE, 3) # checked manually

plot_4 <- ggplot(data = df_info_genome, aes(x = TE_cov_perc, 
                                            y = Model, 
                                            color = Model)) +  
  geom_point(aes(size = 2.5, shape = as.factor(ShapeType))) +
  scale_color_manual(values = map_colorcode) + 
  scale_y_discrete(limits=rev) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 9),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position="none") + 
  labs(x = "TE content (%)", y = "3D model", tag = "D") + 
  annotate(geom="text", x=Inf,y=Inf,
           label=sprintf("italic('p')~'%s'",p_value_TE),
           parse=TRUE, hjust = 1.1, vjust = 1.5, size = 4, color = "black") 
plot_4

# No association between TE content and 3D model 

################################################################################

nested <- (plot_1 | plot_2)/(plot_3|plot_4) + 
  plot_layout(widths = c(1, 1, 1, 1))
nested

ggsave(nested, 
       file=paste0(plotPATH,"Figure3_associations.pdf"), 
       width = 25, height = 15, units = "cm")
