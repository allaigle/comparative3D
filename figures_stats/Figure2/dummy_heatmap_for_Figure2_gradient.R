#!/bin/R

# Author: Alice Laigle

#install.packages("ggplot2")  # Run this line if ggplot2 is not installed
library(ggplot2)

## plot output path (local)
plotPATH <- "/Users/alicelaigle/PhD/fungalHiC/5_final_results/Figure2_genetic_contacts/"

# Sample data creation
set.seed(123)  # For reproducibility
data_matrix <- matrix(rnorm(100), nrow = 10)  # 10x10 matrix of random values
# color palette from Mol* Viewer 
colors <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090",
            "#ffffbf", "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")

# Convert to data frame for ggplot
data_long <- as.data.frame(as.table(data_matrix))
colnames(data_long) <- c("X", "Y", "Value")
# Create the heatmap
g <- ggplot(data_long, aes(x = X, y = Y, fill = Value)) +
  geom_tile() +
  scale_fill_gradientn(colors = colors, name = "Value") + 
  theme_minimal() +
  labs(title = "Heatmap with Color Gradient Legend") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "right"
  )
ggsave(g, 
       file=paste0(plotPATH,"dummy_heatmap_for_Figure2_gradient.pdf"), 
       width = 10, height = 10, units = "cm")