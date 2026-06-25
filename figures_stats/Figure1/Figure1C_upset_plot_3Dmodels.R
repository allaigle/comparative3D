#!/usr/bin/env Rscript

# Goal: Plot panel C - Upset plot of 3D genome architecures and their counts
## Author: Alice Laigle

################################################################################
# Import libraries
library(ggplot2) # v4.0.0
library(tidyverse, warn.conflicts = FALSE) # v2.0.0
library(ggupset) # v0.4.1
theme_set(theme_bw()) # set nice theme

inPATH <- "/Users/alicelaigle/PhD/fungalHiC/0_data/0_collection/"
## plot output path
plotPATH <- "/Users/alicelaigle/PhD/fungalHiC/5_final_results/Figure1_tree_models/"

df <- read.delim(paste0(inPATH,"categorized_models.palette4.txt"),
                 header = TRUE)

df <- df %>%
  mutate(across(c("CT", "Rabl", "Bean"), ~ ifelse(. == 1, TRUE, FALSE)))

combinations_dat <- df  %>%
  mutate(
    combination = pmap(
      list(CT, Rabl, Bean),
      \(lgl1, lgl2, lgl3) {
        c('Chromosome Territories', 'Rabl', 'Bean shape')[c(lgl1, lgl2, lgl3)] 
      }
    )
  )

combinations_dat <- combinations_dat %>%
  mutate(color = case_when(
    combination == "Rabl" ~ "#99CCCC",
    combination == "Bean shape" ~ "#115B75",
    combination == "Chromosome Territories" ~ "#990099",
    TRUE ~ "grey"  # Default color for unmatched combinations
  ))

upset_3D <- combinations_dat %>%
  ggplot(aes(x = combination, fill = color)) +  # Use the new color column
  geom_bar() +
  scale_fill_identity() +  # Use the identity scale for colors
  scale_x_upset()

ggsave(upset_3D, 
       file=paste0(plotPATH,"Figure1_upset_3Dmodels.55_species.pdf"), 
       width = 10, height = 7, units = "cm")

# Note: colors of intermediates have been changed manually, as well has the 
## "Species number" label