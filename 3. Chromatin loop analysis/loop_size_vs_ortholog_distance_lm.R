# Plot distances between coleoid orthologs of loop anchor genes, check R squared for correlation

rm(list = ls())

# Load libraries
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(dplyr)
library(tidyr)

# Reading the files of loop size and orthologous distances
loop_size_dist_eupsc <- read.delim("eupsc_loops_50k+100k.genes_rm_dups_octbi_loop_size_dist.tsv", header = TRUE)
loop_size_dist_sepof <- read.delim("sepof_loops_50k+100k.genes_rm_dups_octbi_loop_size_dist.tsv", header = TRUE)
loop_size_dist_octbi_esc_chr <- read.delim("octbi_loops_50k+100k.genes_rm_dups_eupsc_loop_size_dist.tsv", header = TRUE)
loop_size_dist_octbi_sof_chr <- read.delim("octbi_loops_50k+100k.genes_rm_dups_sepof_loop_size_dist.tsv", header = TRUE)

# Remove rows with NA values
loop_size_dist_eupsc <- na.omit(loop_size_dist_eupsc)
loop_size_dist_sepof <- na.omit(loop_size_dist_sepof)
loop_size_dist_octbi_esc_chr <- na.omit(loop_size_dist_octbi_esc_chr)
loop_size_dist_octbi_sof_chr <- na.omit(loop_size_dist_octbi_sof_chr)

# Remove specific outliers that are missassemblies
loop_size_dist_eupsc <- loop_size_dist_eupsc %>%
  filter(!(loop_size %in% c(14500000, 16700000, 12900000)))

loop_size_dist_octbi_esc_chr <- loop_size_dist_octbi_esc_chr %>%
  filter(!(loop_size %in% c(9400000)))

loop_size_dist_octbi_sof_chr <- loop_size_dist_octbi_sof_chr %>%
  filter(!(loop_size %in% c(9400000)))

# Scatterplot of loop size vs intergenic distance
loop_size_vs_dist_all <- ggplot() +
  geom_point(data = loop_size_dist_eupsc, aes(x = loop_size, y = octbi_average_intergenic_distance, color = "E. scolopes loops, O. bimaculoides orthologs"), size = 3, alpha = 0.7, shape = 16) +
  geom_point(data = loop_size_dist_sepof, aes(x = loop_size, y = octbi_average_intergenic_distance, color = "S. officinalis loops, O. bimaculoides orthologs"), size = 3, alpha = 0.7, shape = 15) +
  geom_point(data = loop_size_dist_octbi_esc_chr, aes(x = loop_size, y = eupsc_average_intergenic_distance, color = "O. bimaculoides loops, E. scolopes orthologs"), size = 3, alpha = 0.7, shape = 17) +
  geom_point(data = loop_size_dist_octbi_sof_chr, aes(x = loop_size, y = sepof_average_intergenic_distance, color = "O. bimaculoides loops, S. officinalis orthologs"), size = 3, alpha = 0.7, shape = 18) +
  labs(
    x = expression(paste("Loop size (bp)")),
    y = expression(paste("Intergenic distance for orthologous gene pairs (bp)"))
  ) +
  scale_color_manual(
    values = c(
      "E. scolopes loops, O. bimaculoides orthologs" = "#117733",
      "S. officinalis loops, O. bimaculoides orthologs" = "#44AA99",
      "O. bimaculoides loops, E. scolopes orthologs" = "#882255",
      "O. bimaculoides loops, S. officinalis orthologs" = "plum3"
    ),
    labels = c(
      "E. scolopes loops, O. bimaculoides orthologs" = expression(italic("E. scolopes") ~ "loops," ~ italic("O. bimaculoides") ~ "orthologs"),
      "S. officinalis loops, O. bimaculoides orthologs" = expression(italic("S. officinalis") ~ "loops," ~ italic("O. bimaculoides") ~ "orthologs"),
      "O. bimaculoides loops, E. scolopes orthologs" = expression(italic("O. bimaculoides") ~ "loops," ~ italic("E. scolopes") ~ "orthologs"),
      "O. bimaculoides loops, S. officinalis orthologs" = expression(italic("O. bimaculoides") ~ "loops," ~ italic("S. officinalis") ~ "orthologs")
    )
  ) +
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.text = element_text(hjust = 0, size = 12),
    legend.title = element_blank()
  )

# Fit linear models
model_eupsc <- lm(octbi_average_intergenic_distance ~ loop_size, data = loop_size_dist_eupsc)
model_sepof <- lm(octbi_average_intergenic_distance ~ loop_size, data = loop_size_dist_sepof)
model_octbi_esc_chr <- lm(eupsc_average_intergenic_distance ~ loop_size, data = loop_size_dist_octbi_esc_chr)
model_octbi_sof_chr <- lm(sepof_average_intergenic_distance ~ loop_size, data = loop_size_dist_octbi_sof_chr)

# Extract R-squared values
r2_eupsc <- summary(model_eupsc)$r.squared
r2_sepof <- summary(model_sepof)$r.squared
r2_octbi_esc_chr <- summary(model_octbi_esc_chr)$r.squared
r2_octbi_sof_chr <- summary(model_octbi_sof_chr)$r.squared

# Annotate R² values in the plot
loop_size_vs_dist_all <- loop_size_vs_dist_all +
  annotate("text", x = Inf, y = Inf, label = paste("R²:", round(r2_eupsc, 2)), hjust = 1.6, vjust = 3.5, size = 5, color = "#117733", fontface = "bold") +
  annotate("text", x = Inf, y = Inf, label = paste("R²:", sprintf("%.2f", r2_sepof)), hjust = 1.6, vjust = 5.5, size = 5, color = "#44AA99", fontface = "bold") +
  annotate("text", x = Inf, y = Inf, label = paste("R²:", round(r2_octbi_esc_chr, 2)), hjust = 1.6, vjust = 7.5, size = 5, color = "#882255", fontface = "bold") +
  annotate("text", x = Inf, y = Inf, label = paste("R²:", round(r2_octbi_sof_chr, 2)), hjust = 1.6, vjust = 9.5, size = 5, color = "plum3", fontface = "bold")

# Print and save the plot
print(loop_size_vs_dist_all)

ggsave("loop_size_vs_dist_all_species.tiff", 
       loop_size_vs_dist_all, dpi = 400, width = 12, height = 7, units = "in")

