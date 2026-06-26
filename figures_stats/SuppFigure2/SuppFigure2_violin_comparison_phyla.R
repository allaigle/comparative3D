#!/usr/bin/env Rscript

# Goal: Plot Supplemental Figure  2
## Violin plot of Asco/Basidiomycota comparisons for genome size, gene/TE contents
## Author: Alice Laigle

################################################################################

### Import libraries
library(ggplot2) # v4.0.0
library(dplyr) # v1.1.4
library(patchwork) # v1.3.2
library(reshape2) # v1.4.4
theme_set(theme_bw()) # set nice theme

################################################################################
# local
plotPATH <- "/Users/alicelaigle/PhD/fungalHiC/5_final_results/SupplementalFigues/"

df_info_genome <- read.delim("/Users/alicelaigle/PhD/fungalHiC/0_data/0_collection/info_genome_TE_SuppFigures2_4.txt",
                        header = TRUE)

################################################################################
### Welch Two Sample t-test

# subset by phylum
df_info_asco <- filter(df_info_genome, Phylum == "Ascomycota")
df_info_basidio <- filter(df_info_genome, Phylum == "Basidiomycota")

# t-test
t_GS <- t.test(df_info_asco$GenomeSize_bp_wo500, df_info_basidio$GenomeSize_bp_wo500) #  0.01298
t.test(df_info_asco$GC_percent_wo500_QUAST, df_info_basidio$GC_percent_wo500_QUAST) # 0.8573
t.test(df_info_asco$Chromosome_Number, df_info_basidio$Chromosome_Number) # 0.8551
t_Gperc <- t.test(df_info_asco$Gene_cov_perc, df_info_basidio$Gene_cov_perc) # 0.002117
t_Tperc <- t.test(df_info_asco$TE_cov_perc, df_info_basidio$TE_cov_perc) # 0.005303

# extract values
t_GS_pval <- round(t_GS$p.value,5)
t_Gperc_pval <- round(t_Gperc$p.value,4)
t_Tperc_pval <- round(t_Tperc$p.value,4)

# extract only Ascomycota and Basidiomycota
df_info_phyla <- filter(df_info_genome, Phylum %in% c("Ascomycota", "Basidiomycota"))

palette_phyla <- c(
  "Ascomycota" = "#3e94d1", 
  "Basidiomycota" = "#cc4b78"
)

GenomeSize_Mbp_wo500 <- (df_info_phyla$GenomeSize_bp_wo500)/1000000

# plot
gs <- ggplot(df_info_phyla, aes(x = Phylum, y = GenomeSize_Mbp_wo500, fill = Phylum)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  geom_point(position = position_jitter(width = 0.2)) + 
  scale_fill_manual(values = palette_phyla) + 
  annotate(geom="text", x=Inf,y=Inf,
           label=sprintf("italic('p')~'%s'",t_GS_pval),
           parse=TRUE, hjust = 1.1, vjust = 1.5, size = 4, color = "black") +
  theme(legend.position = "none") +
  labs(y = "Genome size (Mbp)")
gs

gp <- ggplot(df_info_phyla, aes(x = Phylum, y = Gene_cov_perc, fill = Phylum)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  geom_point(position = position_jitter(width = 0.2)) + 
  scale_fill_manual(values = palette_phyla) + 
  theme(legend.position = "none") +
  annotate(geom="text", x=Inf,y=Inf,
           label=sprintf("italic('p')~'%s'",t_Gperc_pval),
           parse=TRUE, hjust = 1.1, vjust = 1.5, size = 4, color = "black") +
  ylim(0,70) +
  labs(y = "Gene content (%)")
gp

tp <- ggplot(df_info_phyla, aes(x = Phylum, y = TE_cov_perc, fill = Phylum)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  geom_point(position = position_jitter(width = 0.2)) +
  annotate(geom="text", x=Inf,y=Inf,
           label=sprintf("italic('p')~'%s'",t_Tperc_pval),
           parse=TRUE, hjust = 1.1, vjust = 1.5, size = 4, color = "black") +
  scale_fill_manual(values = palette_phyla) +
  ylim(0,100) +
  labs(y = "TE content (%)")
tp

###
nested <- (gs|gp|tp)
nested

ggsave(nested, 
       file=paste0(plotPATH,"SuppFigure2_violin_comparison_phyla.pdf"), 
       width = 30, height = 10, units = "cm")
