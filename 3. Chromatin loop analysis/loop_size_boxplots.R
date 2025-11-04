# Make boxplots of loop sizes

rm(list = ls())

# Load libraries
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(dplyr)
library(tidyr)

# Reading the files
loop_size_eupsc <- read.delim("eupsc_loops_50k+100k.genes_rm_dups_octbi_loop_size_dist.tsv", header = TRUE)
loop_size_sepof <- read.delim("sepof_loops_50k+100k.genes_rm_dups_octbi_loop_size_dist.tsv", header = TRUE)
loop_size_octbi_esc_chr <- read.delim("octbi_loops_50k+100k.genes_rm_dups_eupsc_loop_size_dist.tsv", header = TRUE)
loop_size_octbi_sof_chr <- read.delim("octbi_loops_50k+100k.genes_rm_dups_sepof_loop_size_dist.tsv", header = TRUE)

# Remove specific outliers due to misassemblies
loop_size_eupsc <- loop_size_eupsc %>%
  filter(!(loop_size %in% c(14500000, 16700000, 12900000)))

loop_size_octbi_esc_chr <- loop_size_octbi_esc_chr %>%
  filter(!(loop_size %in% c(9400000)))

loop_size_octbi_sof_chr <- loop_size_octbi_sof_chr %>%
  filter(!(loop_size %in% c(9400000)))

# Filtering for single-ortholog chromosome
loop_size_eupsc_same_chrom <- loop_size_eupsc[loop_size_eupsc$octbi_total_number_of_different_chromosomes_included_in_the_loops_start_and_end == 1, ]
loop_size_sepof_same_chrom <- loop_size_sepof[loop_size_sepof$octbi_total_number_of_different_chromosomes_included_in_the_loops_start_and_end == 1, ]
loop_size_octbi_same_chrom_esc_chr <- loop_size_octbi_esc_chr[loop_size_octbi_esc_chr$eupsc_total_number_of_different_chromosomes_included_in_the_loops_start_and_end == 1, ]
loop_size_octbi_same_chrom_sof_chr <- loop_size_octbi_sof_chr[loop_size_octbi_sof_chr$sepof_total_number_of_different_chromosomes_included_in_the_loops_start_and_end == 1, ]


# Combine data for boxplot
loop_size_combined_same_chrom <- rbind(
  data.frame(loop_size = loop_size_eupsc_same_chrom$loop_size, species = "eupsc"),
  data.frame(loop_size = loop_size_octbi_same_chrom_esc_chr$loop_size, species = "octbi_esc_chr"),
  data.frame(loop_size = loop_size_octbi_same_chrom_sof_chr$loop_size, species = "octbi_sof_chr"),
  data.frame(loop_size = loop_size_sepof_same_chrom$loop_size, species = "sepof")
)

# Reorder factor levels
loop_size_combined_same_chrom$species <- factor(loop_size_combined_same_chrom$species, 
                                                levels = c("eupsc", "sepof", "octbi_esc_chr", "octbi_sof_chr"))

# Boxplot for same-chromosome orthologs
loop_size_box_same_ortho_chrom <- ggplot(loop_size_combined_same_chrom, aes(x = species, y = loop_size, fill = species)) +
  geom_boxplot(
    outlier.shape = 16, 
    outlier.size = 4, 
    alpha = 0.8,
    lwd = 1.2
  ) +
  scale_fill_manual(values = c("eupsc" = "#117733", "sepof" = "#44AA99", "octbi_esc_chr" = "#882255", "octbi_sof_chr" = "plum3")) +
  labs(x = "", y = expression(paste("Loop size (bp)"))) +
  scale_x_discrete(labels = c(
    "eupsc" = expression(italic("E. scolopes")),
    "sepof" = expression(italic("S. officinalis")),
    "octbi_esc_chr" = expression(italic("O. bimaculoides (E. scolopes chr)")),
    "octbi_sof_chr" = expression(italic("O. bimaculoides (S. officinalis chr)"))
  )) +
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.position = "none"
  ) +
  geom_signif(
    comparisons = list(c("eupsc", "sepof"), c("eupsc", "octbi_esc_chr"), c("sepof", "octbi_sof_chr"),c("eupsc", "octbi_sof_chr"), c("sepof", "octbi_esc_chr"), c("octbi_esc_chr", "octbi_sof_chr")),
    map_signif_level = TRUE,
    textsize = 8.5,
    size = 1.5,
    vjust = 0
  )

ggsave("loop_size_box_same_ortho_chrom.tiff", loop_size_box_same_ortho_chrom, dpi = 300, width = 6, height = 7.5, units = "in")

# Repeat the same logic for different-chromosome loops

# Reading the files
loop_size_eupsc_diff_chrom <- loop_size_eupsc[loop_size_eupsc$octbi_total_number_of_different_chromosomes_included_in_the_loops_start_and_end > 1, ]
loop_size_sepof_diff_chrom <- loop_size_sepof[loop_size_sepof$octbi_total_number_of_different_chromosomes_included_in_the_loops_start_and_end > 1, ]
loop_size_octbi_diff_chrom_esc_chr <- loop_size_octbi_esc_chr[loop_size_octbi_esc_chr$eupsc_total_number_of_different_chromosomes_included_in_the_loops_start_and_end > 1, ]
loop_size_octbi_diff_chrom_sof_chr <- loop_size_octbi_sof_chr[loop_size_octbi_sof_chr$sepof_total_number_of_different_chromosomes_included_in_the_loops_start_and_end > 1, ]

loop_size_combined_diff_chrom <- rbind(
  data.frame(loop_size = loop_size_eupsc_diff_chrom$loop_size, species = "eupsc"),
  data.frame(loop_size = loop_size_sepof_diff_chrom$loop_size, species = "sepof"),
  data.frame(loop_size = loop_size_octbi_diff_chrom_esc_chr$loop_size, species = "octbi_esc_chr"),
  data.frame(loop_size = loop_size_octbi_diff_chrom_sof_chr$loop_size, species = "octbi_sof_chr")
)

loop_size_combined_diff_chrom$species <- factor(loop_size_combined_diff_chrom$species, 
                                                levels = c("eupsc", "sepof", "octbi_esc_chr", "octbi_sof_chr"))

# Ensure factor reordering and save the plot
loop_size_box_diff_ortho_chrom <- ggplot(loop_size_combined_diff_chrom, aes(x = species, y = loop_size, fill = species)) +
  geom_boxplot(
    outlier.shape = 16, 
    outlier.size = 4, 
    alpha = 0.8,
    lwd = 1.2
  ) +
  scale_fill_manual(values = c("eupsc" = "#117733", "sepof" = "#44AA99", "octbi_esc_chr" = "#882255", "octbi_sof_chr" = "plum3")) +
  labs(x = "", y = expression(paste("Loop size (bp)")))

ggsave("/Users/users/Desktop/Micro-C/figs_for_paper/Figure4/loop_size_box_diff_ortho_chrom.tiff", loop_size_box_diff_ortho_chrom, dpi = 300, width = 6, height = 7.5, units = "in")

# Repeat the same logic for all loops with genes in, whether they have orthologs or not
# Reading the files
loop_size_eupsc <- read.delim("eupsc_loops_50k+100k.loopsize.tsv.genes_rm_dups", header = TRUE)
loop_size_sepof <- read.delim("sepof_loops_50k+100k.loopsize.tsv.genes_rm_dups", header = TRUE)
loop_size_octbi <- read.delim("octbi_loops_50k+100k.loopsize.tsv.genes_rm_dups", header = TRUE)


# Remove specific outliers due to assemblies
loop_size_eupsc <- loop_size_eupsc %>%
  filter(!(loop_size %in% c(14500000, 16700000, 12900000)))

loop_size_octbi <- loop_size_octbi%>%
  filter(!(loop_size %in% c(9400000)))

loop_size_combined_all_genes <- rbind(
  data.frame(loop_size = loop_size_eupsc$loop_size, species = "eupsc"),
  data.frame(loop_size = loop_size_sepof$loop_size, species = "sepof"),
  data.frame(loop_size = loop_size_octbi$loop_size, species = "octbi")
)

loop_size_combined_all_genes$species <- factor(loop_size_combined_all_genes$species, 
                                                levels = c("eupsc", "sepof", "octbi"))

# Ensure factor reordering and save the plot
loop_size_box_all_genes <- ggplot(loop_size_combined_all_genes, aes(x = species, y = loop_size, fill = species)) +
  geom_boxplot(
    outlier.shape = 16, 
    outlier.size = 4, 
    alpha = 0.8,
    lwd = 1.2
  ) +
  scale_fill_manual(values = c("eupsc" = "#9FE2BF", "sepof" = "#B3EFB2", "octbi" = "#6B5CAE")) +
  labs(x = "", y = expression(paste("Loop size (bp)"))) +
  scale_x_discrete(labels = c(
    "eupsc" = expression(italic("E. scolopes")),
    "sepof" = expression(italic("S. officinalis")),
    "octbi" = expression(italic("O. bimaculoides"))
  )) +
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 18),  
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 16),  
    legend.position = "none"
  ) +
  geom_signif(
    comparisons = list(c("eupsc", "sepof"), c("eupsc", "octbi"), c("sepof", "octbi")),
    map_signif_level = TRUE,
    textsize = 8.5,
    size = 1.5,
    vjust = 0.7,
    y_position = c(21000000, 22000000, 23000000)  # Adjust to prevent overlap
  )

ggsave("loop_size_box_all_with_genes.tiff", loop_size_box_all_genes, dpi = 300, width = 8, height = 8, units = "in")

# Calculate median and mean loop size for each species
summary_stats_same <- loop_size_combined_same_chrom %>%
  group_by(species) %>%
  summarise(
    median_loop_size = median(loop_size, na.rm = TRUE),
    mean_loop_size = mean(loop_size, na.rm = TRUE)
  )

summary_stats_diff <- loop_size_combined_diff_chrom %>%
  group_by(species) %>%
  summarise(
    median_loop_size = median(loop_size, na.rm = TRUE),
    mean_loop_size = mean(loop_size, na.rm = TRUE)
  )

summary_stats <- loop_size_combined_all_genes %>%
  group_by(species) %>%
  summarise(
    median_loop_size = median(loop_size, na.rm = TRUE),
    mean_loop_size = mean(loop_size, na.rm = TRUE)
  )

print(summary_stats_same)
print(summary_stats_diff)
print(summary_stats)


# Assign "Same Chromosome" and "Different Chromosome" labels
loop_size_octbi_same_chrom_esc_chr$type <- "Same Chromosome"
loop_size_octbi_diff_chrom_esc_chr$type <- "Different Chromosome"
loop_size_octbi_same_chrom_sof_chr$type <- "Same Chromosome"
loop_size_octbi_diff_chrom_sof_chr$type <- "Different Chromosome"
loop_size_eupsc_same_chrom$type <- "Same Chromosome"
loop_size_eupsc_diff_chrom$type <- "Different Chromosome"
loop_size_sepof_same_chrom$type <- "Same Chromosome"
loop_size_sepof_diff_chrom$type <- "Different Chromosome"

# Combine data for Octopus bimaculoides (E. scolopes chr)
loop_size_octbi_combined_esc_chr <- rbind(
  data.frame(loop_size = loop_size_octbi_same_chrom_esc_chr$loop_size, type = loop_size_octbi_same_chrom_esc_chr$type, species = "octbi_esc_chr"),
  data.frame(loop_size = loop_size_octbi_diff_chrom_esc_chr$loop_size, type = loop_size_octbi_diff_chrom_esc_chr$type, species = "octbi_esc_chr")
)

# Combine data for Octopus bimaculoides (S. officinalis chr)
loop_size_octbi_combined_sof_chr <- rbind(
  data.frame(loop_size = loop_size_octbi_same_chrom_sof_chr$loop_size, type = loop_size_octbi_same_chrom_sof_chr$type, species = "octbi_sof_chr"),
  data.frame(loop_size = loop_size_octbi_diff_chrom_sof_chr$loop_size, type = loop_size_octbi_diff_chrom_sof_chr$type, species = "octbi_sof_chr")
)

# Combine data for Euprymna scolopes
loop_size_eupsc_combined <- rbind(
  data.frame(loop_size = loop_size_eupsc_same_chrom$loop_size, type = loop_size_eupsc_same_chrom$type, species = "eupsc"),
  data.frame(loop_size = loop_size_eupsc_diff_chrom$loop_size, type = loop_size_eupsc_diff_chrom$type, species = "eupsc")
)

# Combine data for Sepia officinalis
loop_size_sepof_combined <- rbind(
  data.frame(loop_size = loop_size_sepof_same_chrom$loop_size, type = loop_size_sepof_same_chrom$type, species = "sepof"),
  data.frame(loop_size = loop_size_sepof_diff_chrom$loop_size, type = loop_size_sepof_diff_chrom$type, species = "sepof")
)

# Load required library for arranging plots
library(ggpubr)
# Generate individual plots as objects
plot_octbi_esc_chr <- ggplot(loop_size_octbi_combined_esc_chr, aes(x = type, y = loop_size, fill = type)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 4, alpha = 0.8, lwd = 1.2) +
  scale_fill_manual(values = c("Same Chromosome" = "#882255", "Different Chromosome" = "#882255")) +
  labs(x = "", y = expression(paste("Loop size in ", italic("O. bimaculoides"), " (bp)"))) +
  scale_x_discrete(labels = c("Same Chromosome" = "Orthologs on the same E.\nscolopes chromosomes", 
                              "Different Chromosome" = "Orthologs on different E.\nscolopes chromosomes")) +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 18), 
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Rotated x-axis labels
        axis.text.y = element_text(size = 16), 
        legend.position = "none") +
  geom_signif(comparisons = list(c("Same Chromosome", "Different Chromosome")), 
              map_signif_level = TRUE, textsize = 8.5, size = 1.5, vjust = 0, y_position = 21000000)

plot_octbi_sof_chr <- ggplot(loop_size_octbi_combined_sof_chr, aes(x = type, y = loop_size, fill = type)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 4, alpha = 0.8, lwd = 1.2) +
  scale_fill_manual(values = c("Same Chromosome" = "plum3", "Different Chromosome" = "plum3")) +
  labs(x = "", y = expression(paste("Loop size in ", italic("O. bimaculoides"), " (bp)"))) +
  scale_x_discrete(labels = c("Same Chromosome" = "Orthologs on the same S.\nofficinalis chromosomes", 
                              "Different Chromosome" = "Orthologs on different S.\nofficinalis chromosomes")) +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 18), 
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 16), 
        legend.position = "none") +
  geom_signif(comparisons = list(c("Same Chromosome", "Different Chromosome")), 
              map_signif_level = TRUE, textsize = 8.5, size = 1.5, vjust = 0, y_position = 21000000)

plot_eupsc <- ggplot(loop_size_eupsc_combined, aes(x = type, y = loop_size, fill = type)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 4, alpha = 0.8, lwd = 1.2) +
  scale_fill_manual(values = c("Same Chromosome" = "#117733", "Different Chromosome" = "#117733")) +
  labs(x = "", y = expression(paste("Loop size in ", italic("E. scolopes"), " (bp)"))) +
  scale_x_discrete(labels = c("Same Chromosome" = "Orthologs on the same O.\nbimaculoides chromosomes", 
                              "Different Chromosome" = "Orthologs on different O.\nbimaculoides chromosomes")) +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 18), 
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 16), 
        legend.position = "none") +
  geom_signif(comparisons = list(c("Same Chromosome", "Different Chromosome")), 
              map_signif_level = TRUE, textsize = 8.5, size = 1.5, vjust = 0, y_position = 21000000)

plot_sepof <- ggplot(loop_size_sepof_combined, aes(x = type, y = loop_size, fill = type)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 4, alpha = 0.8, lwd = 1.2) +
  scale_fill_manual(values = c("Same Chromosome" = "#44AA99", "Different Chromosome" = "#44AA99")) +
  labs(x = "", y = expression(paste("Loop size in ", italic("S. officinalis"), " (bp)"))) +
  scale_x_discrete(labels = c("Same Chromosome" = "Orthologs on the same O.\nbimaculoides chromosomes", 
                              "Different Chromosome" = "Orthologs on different O.\nbimaculoides chromosomes")) +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 18), 
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 16), 
        legend.position = "none") +
  geom_signif(comparisons = list(c("Same Chromosome", "Different Chromosome")), 
              map_signif_level = TRUE, textsize = 8.5, size = 1.5, vjust = 0, y_position = 21000000)

# Generate and save plots for each species
ggsave("loop_size_octbi_box_same_vs_diff_esc_chr.tiff", 
       plot_octbi_esc_chr, dpi = 300, width = 4, height = 7.2, units = "in")

ggsave("oop_size_octbi_box_same_vs_diff_sof_chr.tiff", 
       plot_octbi_sof_chr, dpi = 300, width = 4, height = 7.2, units = "in")

ggsave("loop_size_eupsc_box_same_vs_diff_oct_chr.tiff", 
       plot_eupsc, dpi = 300, width = 4, height = 7.2, units = "in")

ggsave("loop_size_sepof_box_same_vs_diff_oct_chr.tiff", 
       plot_sepof, dpi = 300, width = 4, height = 7.2, units = "in")

# Arrange the four plots in a 2x2 layout
combined_plot <- ggarrange(
  plot_eupsc, plot_sepof, plot_octbi_esc_chr, plot_octbi_sof_chr, 
  labels = c("A", "B", "C", "D"),
  ncol = 2, nrow = 2
)

# Save the combined plot
ggsave("loop_size_box_same_diff_comparisons_combined.tiff", 
       combined_plot, dpi = 200, width = 10, height = 13, units = "in")


# Print loop sizes
print(summary_stats_same)
print(summary_stats_diff)
print(summary_stats)


