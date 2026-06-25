#!/usr/bin/env Rscript

# Goal: Plot panel 2 - Genetic and interactions descriptions
## Author: Alice Laigle

################################################################################

### Import libraries
library(readr) # v2.1.5
library(ape) # v5.8-1
library(ggplot2) # v4.0.0
library(ggridges) # v0.5.7
library(dplyr) # v1.1.4
library(ggpubr) # v0.6.1
library(patchwork) # v1.3.2
library(reshape2) # v1.4.4 
library(scales) # v1.4.0
library(viridis) # v0.6.5 
library(forcats) # v1.0.1
theme_set(theme_bw()) # set nice theme
source("/Users/alicelaigle/PhD/fungalHiC/3_scripts/variables_names_colors_Fig1-3.R")

################################################################################

### Set PATHS
basePATH <- "/Users/alicelaigle/PhD/fungalHiC" # LOCAL

## path to phylo_tree_order_allDarwin.55.txt, df_info_genome.csv & categorized_models.txt
collectionPATH <- paste0(basePATH,"/0_data/0_collection/")
## path of summaries: TE classes and gene coverage
summAnnoPATH <- paste0(basePATH,"/2_analysis/1_annotations/7_global_summaries/")
## path of summary of types contacts 
contactsPATH <- paste0(basePATH,"/1_dataTreatment/3_conversions/")
## plot output path
plotPATH <- paste0(basePATH,"/5_final_results/Figure2_genetic_contacts/")

################################################################################

### Load files

df_phyloOrder <- read.table(paste0(collectionPATH, "phylo_speciesOrder_55.txt"),
                            sep = "\t", header= TRUE)
df_models <- read.delim(paste0(collectionPATH,"categorized_models.palette4.txt"),
                        header = TRUE)
df_info_genome <- read.table(paste0(collectionPATH,"info_perGenome.55.csv"),
                                sep = ",", header = TRUE)
df_TE_classes <- read.table(paste0(summAnnoPATH,"allDarwin.summary_TE_classes.txt"), 
                            sep = "\t", header = TRUE)
df_interactions <- read.table(paste0(contactsPATH,"allDarwin.summary_interactions.50000.corrected_ICE.txt"), 
                              sep = "\t", header = TRUE)

################################################################################
#### Prepare data 

phyla_order <- c("Basidiomycota", "Ascomycota", "Mucoromycota")
#####################################
### 3D models

subdf <- df_models[,c("Species","CT","Rabl","Bean")]
# Melt df
subdf_models <- melt(subdf, id.vars = "Species")
subdf_models <- subdf_models %>%
  left_join(df_models %>% select(Species, ColorCode), by = "Species") %>%
  mutate(ColorCode = ifelse(value == 0, "white", ColorCode))

subdf_models <- subdf_models %>%
  left_join(df_models %>% select(Species, Full_Species_Name), by = "Species")
subdf_models$Full_Species_Name <- factor(subdf_models$Full_Species_Name, 
                                         levels = df_phyloOrder$Full_Species_Name)

# For main 3D models: bigger circles
color_codes_to_grow_size <- c("#990099", "#99CCCC", "#08085E") 
subdf_models$PointSize <- ifelse(subdf_models$ColorCode %in% color_codes_to_grow_size, 0.5, 0.2)
# For intermediate 3D models: linked small triangles
color_codes_to_connect <- c("#115B75", "#560A7F", "#9292b0")
subdf_models$ShapeType <- ifelse(subdf_models$ColorCode %in% color_codes_to_connect, 17, 16)

#####################################
### info_genome

# Define full name species and phylum as factor 
df_info_genome$Full_Species_Name <- factor(df_info_genome$Full_Species_Name, 
                                           levels = df_phyloOrder$Full_Species_Name)
df_info_genome$Phylum <- as.factor(df_info_genome$Phylum)

#####################################
### TE 
# Add phylum
df_TE_classes <- df_TE_classes %>%
  left_join(select(df_info_genome, Species, Phylum), by = "Species")

# Remove the 5 species out of the 1 million valid contacts' threshold
df_TE_classes <- df_TE_classes %>% 
  filter(!Species %in% dropped_species) # dropped_species is sourced

# Include "Uncertain" in "Unknown"
df_TE_classes <- df_TE_classes %>%
  mutate(Category = recode(Category, 
                           "Uncertain" = "Unknown"))
df_TE_classes <- df_TE_classes %>%
  mutate(Category = recode(Category, 
                           "Retroposon" = "LTR")) # really few RT

# Create palette of TE types
palClasses <- c("DNA" = "#000000", "LINE" = "#E69F00", "LTR" = "#56B4E9", 
                "PLE" = "#009E73", "RC" = "#F0E442", "
                SINE" = "#D55E00", "Unknown" = "#999999")

# define as factor
df_TE_classes$Category <- as.factor(df_TE_classes$Category)
df_TE_classes$Species <- factor(df_TE_classes$Species, levels = df_phyloOrder$Species)
df_TE_classes$Phylum <- as.factor(df_TE_classes$Phylum)

#####################################
### Interactions
# Subset to reorder
mdf_interactions <- df_interactions[,c("Species",
                                         "short_cis_interactions",
                                         "long_cis_interactions",
                                         "trans_interactions")] # subset
# Melt df
mdf_interactions <- melt(mdf_interactions, id.vars = "Species")
mdf_interactions$Species <- factor(mdf_interactions$Species, 
                                     levels = df_phyloOrder$Species)
# reorder interactions
mdf_interactions <- mdf_interactions %>%
  mutate(variable = fct_relevel(variable, "trans_interactions" , 
                                "long_cis_interactions", 
                                "short_cis_interactions"))
# Add phylum
mdf_interactions <- mdf_interactions %>%
  left_join(select(df_info_genome, Species, Phylum), by = "Species")

# Define as factor
mdf_interactions$Species <- factor(mdf_interactions$Species, 
                              levels = df_phyloOrder$Species) 
mdf_interactions$Phylum <- as.factor(mdf_interactions$Phylum)

################################################################################

### 1. Dot plot of 3D models as in Figure 1

g_models <- ggplot(subdf_models, 
                   aes(x = Full_Species_Name, y = variable, color = ColorCode)) + 
  geom_point(aes(size = PointSize, shape = as.factor(ShapeType)), stroke = 0.1) + 
  geom_line(data = subdf_models %>% filter(ColorCode %in% color_codes_to_connect), 
            aes(group = Species), linewidth = 0.5) +
  scale_color_identity() +
  scale_size_continuous(range = c(2, 4.5)) +
  scale_shape_manual(values = c(16, 17)) +
  theme(legend.position = "none", 
        axis.title = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 0.2, color = "black"),
        panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size=9, vjust = 0.5),
        axis.text.y = element_text(face = "italic", size = 10)) +
  coord_flip() +
  labs(x = "Species", y = "3D model")
g_models

################################################################################
### 2. Dot plot of genome size filled with GC content as gradient

g_GSize <- ggplot(df_info_genome, aes(x = Full_Species_Name, y = GenomeSize_bp_wo500)) + 
  geom_point(aes(fill = GC_percent_wo500_QUAST), 
             size = 5, shape = 21, color = "black", stroke = 0.3) +
  scale_fill_gradient(low = "#E7DDDA", high = "#6D3625", name = "GC %") +
  theme(axis.title = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(colour = "black", size=9, vjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(size = 0.2, color = "black"),
        legend.position = "bottom") + 
  coord_flip() + 
  scale_y_log10(breaks = c(2e+7, 5e+7, 1e+8)) +
  labs(y = "Genome size (bp)")
g_GSize

(min(df_info_genome$GenomeSize_bp_wo500))/1000000 # in Mb
(max(df_info_genome$GenomeSize_bp_wo500))/1000000
(median(df_info_genome$GenomeSize_bp_wo500))/1000000

min(df_info_genome$GC_percent_wo500_QUAST)
max(df_info_genome$GC_percent_wo500_QUAST)
median(df_info_genome$GC_percent_wo500_QUAST)

################################################################################

### 3. Dot plot of chromosome number

g_chrNumber <- ggplot(df_info_genome, aes(x=Full_Species_Name, y=Chromosome_Number)) + 
  geom_point(fill="#125C3B",size=4, shape = 21, color = "black", 
             stroke = 0.3) +
  coord_flip() + 
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 9),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position="bottom",
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(size = 0.2, color = "black"),
        legend.box.background = element_rect(colour = "black")) +
  labs(y="Chromosome number") 

g_chrNumber

min(df_info_genome$Chromosome_Number)
max(df_info_genome$Chromosome_Number)
median(df_info_genome$Chromosome_Number)

################################################################################

### 4. Bar plot of TE content

gTEclasses <- ggplot(df_TE_classes, 
                     aes(fill=Category, y=Coverage_percent, x=Species)) + 
  geom_bar(aes(fill = Category), position="stack", stat="identity") +
  scale_fill_manual(values = palClasses) + 
  theme(axis.text.x= element_text(colour = "black", 
                                  size = 9,
                                  vjust = 0.5),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(size = 0.2, color = "black"),
        legend.position = "bottom") +
  coord_flip() + 
  labs(y= "TE coverage (%)", fill="Category")

gTEclasses

min(df_info_genome$TE_cov_perc)
max(df_info_genome$TE_cov_perc)
median(df_info_genome$TE_cov_perc)

################################################################################

### 5. Stacked barplot of types of contacts 

g_contacts <- ggplot(mdf_interactions, aes(x = Species, y=value, fill = variable)) +
  geom_bar(position = "fill", stat = "identity", color = "black", 
           linewidth = 0.1) +
  scale_fill_manual(values = c("short_cis_interactions" = "#8a5e50", 
                               "long_cis_interactions" = "#E7DDDA", 
                               "trans_interactions" = "#709d89"),
                    labels = c("long_cis_interactions" = expression("long " * italic("cis")),
                               "short_cis_interactions" = expression("short " * italic("cis")),
                               "trans_interactions" = expression(italic("trans")))) +
  scale_y_continuous(labels = label_percent(accuracy = 1, scale = 100, suffix = "")) + 
  labs(x = "Species ordered by Phylogenetic", y = "Valid Hi-C interactions (%)", 
       fill = "Type") +
  theme(axis.title = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x= element_text(colour = "black", 
                                  size = 9,
                                  angle = 0, 
                                  vjust = 0.5),
        axis.line = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 0.2, color = "black"),
        panel.background = element_blank(),
        legend.position="bottom") +
  coord_flip() 
g_contacts

min(df_interactions$cis_perc)
max(df_interactions$cis_perc)
median(df_interactions$cis_perc)
min(df_interactions$trans_perc)
max(df_interactions$trans_perc)
median(df_interactions$trans_perc)

################################################################################

### Plot them all at once

nested <- (g_models | g_GSize | g_chrNumber | gTEclasses | g_contacts) + 
  plot_layout(widths = c(0.5, 0.75, 0.5, 1.5, 1))

nested

# save
ggsave(nested, 
       file=paste0(plotPATH,"Figure2_general_description.55_species.v2.pdf"), 
       width = 25, height = 25, units = "cm")