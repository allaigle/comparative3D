#!/bin/R

# Goal: Plot the Supplemental Figure 4 - phylogenetic regressions and ANOVAs

# Author: Alice Laigle
# Adapted from Camille Cornet's script

library(ape) # for function corPagel - v5.8-1
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

## Add columns from genome info
df_interactions <- merge(df_interactions, 
                         df_info_genome[, c("Species", "GenomeSize_bp_wo500", 
                                            "TE_cov_perc", "Model", "ColorCode")], 
                         by = "Species", 
                         all.x = TRUE) 

# Calculate cis/trans ratio
df_info_genome$Ratio <- (df_info_genome$Cis_Interactions_perc/df_info_genome$Trans_Interactions_perc)

# Make intermediates triangles
color_codes_to_reshape <- c("#115B75", "#560A7F", "#9292b0")
df_info_genome$ShapeType <- ifelse(df_info_genome$ColorCode %in% color_codes_to_reshape, 17, 16)
df_interactions$ShapeType <- ifelse(df_interactions$ColorCode %in% color_codes_to_reshape, 17, 16)
###############################################################################

# Define a variance-covariance structure
pagel <- corPagel(1, phy = phylo_tree, fixed = FALSE, form = ~Species)

# create function to run model and extract p-value, AIC and BIC
extract_model_stats <- function(data, response, predictor, correlation_structure) {
  formula <- as.formula(paste(response, "~", predictor)) # 
  model <- gls(formula, data = data, correlation = correlation_structure) # run model
  
  p_value <- summary(model)$tTable[2, "p-value"] # extract p-value
  
  model_aic <- AIC(model) # extract AIC
  model_bic <- BIC(model) # extract BIC
  
  results <- list( # Create a list to return results
    p_value = p_value,
    aic = model_aic,
    bic = model_bic
  )
  return(results)
}

################################################################################
# 0. Quality control
## check if the genome size is affecting the number of valid contacts 

res_QC <- extract_model_stats(df_interactions, "total_interactions", 
                    "GenomeSize_bp_wo500", pagel)
res_QC[[1]] # check p-value to adapt the number of numbers rounded - 0.3644896
pval_QC <- round(res_QC[[1]], 3)
res_QC[[2]] # AIC
res_QC[[3]] # BIC 

df_interactions$GenomeSize_Mbp <- (df_interactions$GenomeSize_bp_wo500)/1000000
df_interactions$total_interactions_M <- (df_interactions$total_interactions)/1000000

plot_0 <- ggplot(data = df_interactions, aes(x = GenomeSize_Mbp, 
                                             y = total_interactions_M,
                                             colour = Model)) +  
  geom_point(aes(size = 0.25, shape = as.factor(ShapeType))) +
  scale_color_manual(values = map_colorcode) +
  theme(axis.title = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text = element_text(size = 9),
        legend.position = "none"
  ) + 
  labs(x = "Genome size (Mbp)", y = "Total interactions (Million)", color = "3D model",
       tag = "A") + 
  annotate(geom="text", x = Inf, y = Inf,
           label=sprintf("italic('p')~'%s'",pval_QC),
           parse=TRUE, hjust = 1.1, vjust = 1.5, size = 4, color = "black")
plot_0

# No correlation between the genome size and the number of interactions

################################################################################

## 1. Is the genome size correlated to the chromosome number ?

res_gs_nb <- extract_model_stats(df_info_genome, "GenomeSize_bp_wo500", 
                              "Chromosome_Number", pagel)
res_gs_nb[[1]] # 0.03024583
pval_gs_nb <- round(res_gs_nb[[1]], 4)
res_gs_nb[[2]] # AIC
res_gs_nb[[3]] # BIC 

df_info_genome$GenomeSize_Mbp <- (df_info_genome$GenomeSize_bp_wo500)/1000000

plot_1 <- ggplot(data = df_info_genome, aes(x = GenomeSize_Mbp, 
                                            y = Chromosome_Number,
                                            colour = Model)) +  
  geom_point(aes(size = 0.05, shape = as.factor(ShapeType))) +
  scale_color_manual(values = map_colorcode) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 9),
        legend.position = "none"
  ) + 
  labs(x = "Genome size (Mbp)", y = "Chromosome number", color = "3D model",
       tag = "B") + 
  annotate(geom="text", x = Inf, y = Inf,
           label=sprintf("italic('p')~'%s'",pval_gs_nb),
           parse=TRUE, hjust = 1.1, vjust = 1.5, size = 4, color = "black")
plot_1

# Correlation between the genome size and the chromosome number

################################################################################

## 2. Is the genome size correlated to the gene coverage ?

res_gs_gcov <- extract_model_stats(df_info_genome, "GenomeSize_bp_wo500", 
                                 "Gene_cov_perc", pagel)
res_gs_gcov[[1]] # 1.928756e-13
pval_gs_gcov <- round(res_gs_gcov[[1]], 14)
res_gs_gcov[[2]] # AIC
res_gs_gcov[[3]] # BIC 

plot_2 <- ggplot(data = df_info_genome, aes(x = GenomeSize_Mbp, 
                                            y = Gene_cov_perc,
                                            colour = Model)) +  
  geom_point(aes(size = 0.05, shape = as.factor(ShapeType))) +
  scale_color_manual(values = map_colorcode) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 9),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none"
  ) + 
  labs(x = "Genome size (Mbp)", y = "Gene content (%)", color = "3D model",
       tag = "C") + 
  annotate(geom="text", x = Inf, y = Inf,
           label=sprintf("italic('p')~'%s'",pval_gs_gcov),
           parse=TRUE, hjust = 1.1, vjust = 1.5, size = 4, color = "black")
plot_2

# Very strong correlation between the genome size and the gene coverage

################################################################################

## 3. Is the genome size correlated with CG content ?

res_gs_gc <- extract_model_stats(df_info_genome, "GenomeSize_bp_wo500", 
                                 "GC_percent_wo500_QUAST", pagel)
res_gs_gc[[1]] # 3.428599e-05
pval_gs_gc <- round(res_gs_gc[[1]], 7)
pval_gs_gc
res_gs_gc[[2]] # AIC
res_gs_gc[[3]] # BIC 
plot_3 <- ggplot(data = df_info_genome, aes(x = GenomeSize_Mbp, 
                                            y = GC_percent_wo500_QUAST,
                                            colour = Model)) +  
  geom_point(aes(size = 0.05, shape = as.factor(ShapeType))) +
  scale_color_manual(values = map_colorcode) +
  ylim(25,65) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 9),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none"
  ) + 
  labs(x = "Genome size (Mbp)", y = "GC content (%)", color = "3D model",
       tag = "D") + 
  annotate(geom="text", x = Inf, y = Inf,
           label=sprintf("italic('p')~'%s'",pval_gs_gc),
           parse=TRUE, hjust = 1.1, vjust = 1.5, size = 4, color = "black")
plot_3

# Strong correlation between the genome size and the GC content

################################################################################
## 4. Is the CG content correlated with TE coverage ?

res_te_gc <- extract_model_stats(df_info_genome, "TE_cov_perc", 
                                   "GC_percent_wo500_QUAST", pagel)
res_te_gc[[1]] # 0.001664954
pval_te_gc <- round(res_te_gc[[1]], 4)
pval_te_gc
res_te_gc[[2]] # AIC
res_te_gc[[3]] # BIC 

plot_4 <- ggplot(data = df_info_genome, aes(x = TE_cov_perc, 
                                            y = GC_percent_wo500_QUAST,
                                            colour = Model)) +  
  geom_point(aes(size = 1, shape = as.factor(ShapeType))) +
  scale_color_manual(values = map_colorcode) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 9),
        legend.position = "none"
  ) + 
  labs(x = "GC content (%)", y = "TE content (%)", color = "3D model",
       tag = "E") + 
  annotate(geom="text", x = Inf, y = Inf,
           label=sprintf("italic('p')~'%s'",pval_te_gc),
           parse=TRUE, hjust = 1.1, vjust = 1.5, size = 4, color = "black")
plot_4

# Strong correlation between the GC content and TE coverage

################################################################################

## 5. Is the genome size correlated to specific TEs?

### LTR ########################################################################

res_gs_LTR <- extract_model_stats(df_info_genome, "GenomeSize_bp_wo500", 
                                  "LTR", pagel)
res_gs_LTR[[1]] # 3.233179e-05
pval_gs_LTR <- round(res_gs_LTR[[1]], 7)
pval_gs_LTR
res_gs_LTR[[2]] # AIC
res_gs_LTR[[3]] # BIC 

plot_5 <- ggplot(data = df_info_genome, aes(x = GenomeSize_Mbp, 
                                            y = LTR,
                                            colour = Model)) +  
  geom_point(aes(size = 0.05, shape = as.factor(ShapeType))) +
  scale_color_manual(values = map_colorcode) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 9),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none"
  ) + 
  labs(x = "Genome size (Mbp)", y = "LTR content (%)", color = "3D model",
       tag = "F") + 
  annotate(geom="text", x = Inf, y = Inf,
           label=sprintf("italic('p')~'%s'",pval_gs_LTR),
           parse=TRUE, hjust = 1.1, vjust = 1.5, size = 4, color = "black")
plot_5
# Very strong correlation between the genome size and the LTR category

### DNA ########################################################################

res_gs_DNA <- extract_model_stats(df_info_genome, "GenomeSize_bp_wo500", 
                                 "DNA", pagel)
res_gs_DNA[[1]] # 8.454129e-09
pval_gs_DNA <- round(res_gs_DNA[[1]], 11)
pval_gs_DNA
res_gs_DNA[[2]] # AIC
res_gs_DNA[[3]] # BIC 

plot_6 <- ggplot(data = df_info_genome, aes(x = GenomeSize_Mbp, 
                                            y = DNA,
                                            colour = Model)) +  
  geom_point(aes(size = 0.05, shape = as.factor(ShapeType))) +
  scale_color_manual(values = map_colorcode) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 9),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none"
  ) + 
  labs(x = "Genome size (Mbp)", y = "DNA content (%)", color = "3D model",
       tag = "G") + 
  annotate(geom="text", x = Inf, y = Inf,
           label=sprintf("italic('p')~'%s'",pval_gs_DNA),
           parse=TRUE, hjust = 1.1, vjust = 1.5, size = 4, color = "black")
plot_6
# Very strong correlation between the genome size and the DNA category

### LINE ########################################################################
res_gs_LINE <- extract_model_stats(df_info_genome, "GenomeSize_bp_wo500", 
                                  "LINE", pagel)
res_gs_LINE[[1]] # 0.100565
pval_gs_LINE <- round(res_gs_LINE[[1]], 3)
pval_gs_LINE
res_gs_LINE[[2]] # AIC
res_gs_LINE[[3]] # BIC 

plot_7 <- ggplot(data = df_info_genome, aes(x = GenomeSize_Mbp, 
                                            y = LINE,
                                            colour = Model)) +  
  geom_point(aes(size = 0.05, shape = as.factor(ShapeType))) +
  scale_color_manual(values = map_colorcode) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 9),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none"
  ) + 
  labs(x = "Genome size (Mbp)", y = "LINE content (%)", color = "3D model",
       tag = "H") + 
  annotate(geom="text", x = Inf, y = Inf,
           label=sprintf("italic('p')~'%s'",pval_gs_LINE),
           parse=TRUE, hjust = 1.1, vjust = 1.5, size = 4, color = "black")
plot_7

# No correlation between the genome size and LINE


### PLE ########################################################################
res_gs_PLE <- extract_model_stats(df_info_genome, "GenomeSize_bp_wo500", 
                                   "PLE", pagel)
res_gs_PLE[[1]] # 0.756261
pval_gs_PLE <- round(res_gs_PLE[[1]], 3)
pval_gs_PLE
res_gs_PLE[[2]] # AIC
res_gs_PLE[[3]] # BIC 

plot_8 <- ggplot(data = df_info_genome, aes(x = GenomeSize_Mbp, 
                                            y = PLE,
                                            colour = Model)) +  
  geom_point(aes(size = 0.05, shape = as.factor(ShapeType))) +
  scale_color_manual(values = map_colorcode) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 9),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none"
  ) + 
  labs(x = "Genome size (Mbp)", y = "PLE content (%)", color = "3D model",
       tag = "I") + 
  annotate(geom="text", x = Inf, y = Inf,
           label=sprintf("italic('p')~'%s'",pval_gs_PLE),
           parse=TRUE, hjust = 1.1, vjust = 1.5, size = 4, color = "black")
plot_8
# No correlation between the genome size and PLE


### RC #########################################################################
res_gs_RC <- extract_model_stats(df_info_genome, "GenomeSize_bp_wo500", 
                                  "RC", pagel)
res_gs_RC[[1]] # 0.1009124
pval_gs_RC <- round(res_gs_RC[[1]], 3)
pval_gs_RC
res_gs_RC[[2]] # AIC
res_gs_RC[[3]] # BIC 


plot_9 <- ggplot(data = df_info_genome, aes(x = GenomeSize_Mbp, 
                                            y = RC,
                                            colour = Model)) +  
  geom_point(aes(size = 0.05, shape = as.factor(ShapeType))) +
  scale_color_manual(values = map_colorcode) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 9),
        legend.position = "none"
  ) + 
  labs(x = "Genome size (Mbp)", y = "RC content (%)", color = "3D model",
       tag = "J") + 
  annotate(geom="text", x = Inf, y = Inf,
           label=sprintf("italic('p')~'%s'",pval_gs_RC),
           parse=TRUE, hjust = 1.1, vjust = 1.5, size = 4, color = "black")
plot_9
# No correlation between the genome size and RC


### SINE #######################################################################
res_gs_SINE <- extract_model_stats(df_info_genome, "GenomeSize_bp_wo500", 
                                   "SINE", pagel)
res_gs_SINE[[1]] # 6.253369e-06
pval_gs_SINE <- round(res_gs_SINE[[1]], 8)
pval_gs_SINE
res_gs_SINE[[2]] # AIC
res_gs_SINE[[3]] # BIC 

plot_10 <- ggplot(data = df_info_genome, aes(x = GenomeSize_Mbp, 
                                            y = SINE,
                                            colour = Model)) +  
  geom_point(aes(size = 0.05, shape = as.factor(ShapeType))) +
  scale_color_manual(values = map_colorcode) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 9),
        legend.position = "none"
  ) + 
  labs(x = "Genome size (Mbp)", y = "SINE content (%)", color = "3D model",
       tag = "K") + 
  annotate(geom="text", x = Inf, y = Inf,
           label=sprintf("italic('p')~'%s'",pval_gs_SINE),
           parse=TRUE, hjust = 1.1, vjust = 1.5, size = 4, color = "black")
plot_10

min(df_info_genome$SINE)
max(df_info_genome$SINE)
median(df_info_genome$SINE)

################################################################################

# 11. Is the genome size correlated to the ratio cis/trans?
res_gs_ratio <- extract_model_stats(df_info_genome, "Ratio", 
                                   "GenomeSize_bp_wo500", pagel)
res_gs_ratio[[1]] # 0.03767548
pval_gs_ratio <- round(res_gs_ratio[[1]], 4)
pval_gs_ratio
res_gs_ratio[[2]] # AIC
res_gs_ratio[[3]] # BIC 

plot_11 <- ggplot(data = df_info_genome, aes(x = GenomeSize_Mbp, 
                                             y = Ratio,
                                             colour = Model)) +  
  geom_point(aes(size = 0.05, shape = as.factor(ShapeType))) +
  scale_color_manual(values = map_colorcode) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 9),
        legend.position = "none"
  ) + 
  labs(x = "Genome size (Mbp)", y = "Ratio of cis/trans interactions", 
       color = "3D model", tag = "L") + 
  annotate(geom="text", x = Inf, y = Inf,
           label=sprintf("italic('p')~'%s'",pval_gs_ratio),
           parse=TRUE, hjust = 1.1, vjust = 1.5, size = 4, color = "black")
plot_11

################################################################################

# 11. Is the genome size correlated to the ratio cis/trans?
res_nb_ratio <- extract_model_stats(df_info_genome, "Ratio", 
                                    "Chromosome_Number", pagel)
res_nb_ratio[[1]] # 0.9139268
pval_nb_ratio <- round(res_nb_ratio[[1]], 3)

plot_12 <- ggplot(data = df_info_genome, aes(x = Chromosome_Number, 
                                             y = Ratio,
                                             colour = Model)) +  
  geom_point(aes(size = 0.05, shape = as.factor(ShapeType))) +
  scale_color_manual(values = map_colorcode) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 9),
        legend.position = "none"
  ) + 
  labs(x = "Chromosome Number", y = "Ratio of cis/trans interactions", 
       color = "3D model", tag = "M") + 
  annotate(geom="text", x = Inf, y = Inf,
           label=sprintf("italic('p')~'%s'",pval_nb_ratio),
           parse=TRUE, hjust = 1.1, vjust = 1.5, size = 4, color = "black")
plot_12

################################################################################
#plot_blank <- ggplot() + theme_void() # create blanck plot
nested <- (plot_0 | plot_1 | plot_2) /
  (plot_3 | plot_4 | plot_5) /
  (plot_6 | plot_7 | plot_8) /
  (plot_9 | plot_10 | plot_11) +
  plot_layout(guides = "collect") & theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
nested

ggsave(nested, 
       file=paste0(plotPATH,"SuppFigure4_phyloreg.pdf"), 
       width = 20, height = 25, units = "cm")
