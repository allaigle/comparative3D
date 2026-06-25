#!/usr/bin/env Rscript

# Goal: Plot panel 1 - Dot plot of categorized 3D models without displaying 
#                      species names
## Author: Alice Laigle

################################################################################
# Import libraries
library(reshape2) # v1.4.4
library(ggplot2) # v4.0.0
theme_set(theme_bw()) # set nice theme
source("/Users/alicelaigle/PhD/fungalHiC/3_scripts/variables_names_colors_Fig1-3.R")

# Paths
inPATH <- "/Users/alicelaigle/PhD/fungalHiC/0_data/0_collection/"
## path where the order of species is written (phylo_tree_order_allDarwin.55.txt)
phyloOrderPATH <- "/Users/alicelaigle/PhD/fungalHiC/0_data/0_collection/"
## plot output path
plotPATH <- "/Users/alicelaigle/PhD/fungalHiC/5_final_results/Figure1_tree_models/"

df <- read.delim(paste0(inPATH,"categorized_models.palette4.txt"),
                 header = TRUE)
phyloOrder <- read.table(paste0(phyloOrderPATH, "phylo_speciesOrder_55.txt"),
                         sep = "\t", header= TRUE)

subdf <- df[,c("Species","CT","Rabl","Bean")]

# Melt df
subdf_models <- melt(subdf, id.vars = "Species")
subdf_models <- subdf_models %>%
  left_join(df %>% select(Species, ColorCode), by = "Species") %>%
  mutate(ColorCode = ifelse(value == 0, "white", ColorCode))

subdf_models <- subdf_models %>%
  left_join(df %>% select(Species, Full_Species_Name), by = "Species") 
subdf_models$Full_Species_Name <- factor(subdf_models$Full_Species_Name, 
                                         levels = phyloOrder$Full_Species_Name)

color_codes_to_connect <- c("#115B75", "#560A7F", "#9292b0")
color_codes_to_grow_size <- c("#990099", "#99CCCC", "#08085E")

subdf_models$ShapeType <- ifelse(subdf_models$ColorCode %in% 
                                   color_codes_to_connect, 17, 16)
subdf_models$PointSize <- ifelse(subdf_models$ColorCode %in% 
                                   color_codes_to_grow_size, 3.5, 3)

# Create the ggplot
g_models <- ggplot(subdf_models, 
                   aes(x = Full_Species_Name, y = variable, color = ColorCode)) + 
  geom_point(aes(size = PointSize, shape = as.factor(ShapeType)), stroke = 0.1) + 
  geom_line(data = subdf_models %>% filter(ColorCode %in% color_codes_to_connect), 
            aes(group = Species), linewidth = 0.5) +
  scale_color_identity() +
  scale_size_continuous(range = c(3, 7)) +
  scale_shape_manual(values = c(16, 17)) +
  theme(legend.position = "none", 
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank()) +
  coord_flip() +
  labs(x = "Species", y = "3D model")
g_models

# save
ggsave(g_models, 
       file=paste0(plotPATH,"Figure1_dotPlot_models.55_species.woSpecies.palette4.pdf"), 
       width = 4, height = 35, units = "cm")