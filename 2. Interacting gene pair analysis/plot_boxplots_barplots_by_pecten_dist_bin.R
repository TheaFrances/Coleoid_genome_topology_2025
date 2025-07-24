# Split gene pairs into three groups, <5 Mb, 5–15 Mb, and ≥15 Mb distance between orthologs on P. maximus chromosomes
# Create boxplots of genomic distance between orthologous gene pairs across interaction statuses, coloured by these P. maximus distance groups
# Test for significant differences between genomic distances across *P. maximus* distance bins using BH-corrected Wilcoxon tests
# Then plot a barplot for the number of gene pairs each interaction category and P. maximus distance category

# Clear workspace
rm(list = ls())

# Load required libraries
library(dplyr)
library(ggplot2)
library(tidyr)

# Load data
pec_dists <- read.delim("409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof_PEC_dist.txt")

# Remove NAs (things on different chromosomes)
pec_dists_na_rm <- na.omit(pec_dists)

# Binning based on Pecten_dist
pec_dists_na_rm <- pec_dists_na_rm %>%
  mutate(Bin_Label_3 = case_when(
    Pecten_dist < 5000000 ~ "bin1",
    Pecten_dist >= 5000000 & Pecten_dist < 15000000 ~ "bin2",
    Pecten_dist >= 15000000 ~ "bin3",
    TRUE ~ NA_character_
  ))


# Colors for bins
bin_colors_3 <- c( "yellow3", "cyan4","darkslateblue")

pec_dists_ave <- pec_dists_na_rm

# Save merged file for future analyses----
write.table(pec_dists_ave, "409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_pec_bins_with_sof.txt", append = FALSE, sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE, quote = FALSE)


# Interaction subsets
interaction_types <- c("interacting_all_species", "not_interacting_any_species", "interacting_deca_only", "interacting_octbi_only")
interaction_data <- lapply(interaction_types, function(type) {
  subset(pec_dists_ave, Interaction_status == type)
})

names(interaction_data) <- interaction_types

# Set the order of Interaction_status
pec_dists_ave$Interaction_status <- factor(pec_dists_ave$Interaction_status, 
                                           levels = interaction_types)

# Print the number of genes in each bin category
gene_counts <- pec_dists_ave %>%
  group_by(Bin_Label_3) %>%
  summarize(count = n())
print(gene_counts)

# Subset to make different plots, because they don't all fit in the same one
pec_dists_ave_not_int <- subset(pec_dists_ave, subset=(Interaction_status=="not_interacting_any_species" | Interaction_status=="interacting_octbi_only" ))
pec_dists_ave_int <- subset(pec_dists_ave, subset=(Interaction_status=="interacting_all_species" | Interaction_status=="interacting_deca_only")) 

 # Plotting the boxplots
grouped_box_plot_not_int <- ggplot(pec_dists_ave_not_int, aes(x = Interaction_status, y = Genomic_distance_eupsc_bp, fill = Bin_Label_3)) +
  geom_boxplot(outlier.size = 4, size = 2) +
  labs(x = NULL, y = expression(paste("Genomic distance in ", italic("E. scolopes"), " (bp)"))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 50), # Adjust font size for y-axis labels
        axis.title.y = element_text(size = 50), # Adjust font size for y-axis title
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_fill_manual(name = NULL, values = bin_colors_3, labels = c(
    expression(paste("< 5Mb apart in ", italic("P. maximus"), ", (n = 4464)")),
    expression(paste("≥ 5Mb and < 15Mb apart in ", italic("P. maximus"), ", (n =  6588)")),
    expression(paste("≥ 15Mb apart in ", italic("P. maximus"), ", (n = 7935)"))
  )) + coord_cartesian(y = c(0,200000000))

# Save the plot
ggsave("grouped_box_>=10_pec_bins_ave_eupsc_dist_not_int_with_sof.tiff", 
       grouped_box_plot_not_int, 
       width = 15, 
       height = 30, 
       units = "in", 
       dpi = 400)

# Plotting the boxplots
grouped_box_plot_int <- ggplot(pec_dists_ave_int, aes(x = Interaction_status, y = Genomic_distance_eupsc_bp, fill = Bin_Label_3)) +
  geom_boxplot(outlier.size = 4, size = 2) +
  labs(x = NULL, y = expression(paste("Genomic distance in ", italic("E. scolopes"), " (bp)"))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 50), # Adjust font size for y-axis labels
        axis.title.y = element_text(size = 50), # Adjust font size for y-axis title
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_fill_manual(name = NULL, values = bin_colors_3, labels = c(
    expression(paste("< 5Mb apart in ", italic("P. maximus"), ", (n = 4464)")),
    expression(paste("≥ 5Mb and < 15Mb apart in ", italic("P. maximus"), ", (n =  6588)")),
    expression(paste("≥ 15Mb apart in ", italic("P. maximus"), ", (n = 7935))"))
  )) + coord_cartesian(y = c(0,7500000))

# Save the plot
ggsave("grouped_box_>=10_pec_bins_ave_eupsc_dist_int_with_sof.tiff", 
       grouped_box_plot_int, 
       width = 15, 
       height = 30, 
       units = "in", 
       dpi = 400)

# Perform pairwise Wilcoxon tests within each Interaction_status group
pairwise_results <- lapply(interaction_types, function(type) {
  subset_data <- pec_dists_ave %>% filter(Interaction_status == type)
  
  # Perform the pairwise Wilcoxon test
  test_result <- pairwise.wilcox.test(subset_data$Genomic_distance_eupsc_bp, subset_data$Bin_Label_3, p.adjust.method = "BH")
  
  # Extract p-values and format results
  result_df <- as.data.frame(test_result$p.value)
  
  # Convert the matrix to a long format data frame
  result_long <- result_df %>%
    tibble::rownames_to_column(var = "Group1") %>%
    pivot_longer(cols = -Group1, names_to = "Group2", values_to = "P_Value") %>%
    mutate(Interaction_status = type)
  
  result_long
})

# Combine results into a single data frame
pairwise_results_df <- bind_rows(pairwise_results)

write.table(pairwise_results_df, "pec_bins_eupsc_dist_wilcox_with_sof.txt", append = FALSE, sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Calculate medians for each group
median_results <- pec_dists_ave %>%
  group_by(Interaction_status, Bin_Label_3) %>%
  summarize(Median_Genomic_Distance = median(Genomic_distance_eupsc_bp), .groups = 'drop')

# Print the medians
print(median_results)

# Calculate medians for each group
mean_results <- pec_dists_ave%>%
  group_by(Interaction_status, Bin_Label_3) %>%
  summarize(Mean_Genomic_Distance = mean(Genomic_distance_eupsc_bp), .groups = 'drop')

# Print the medians
print(mean_results)

#Grouped barplot----
pec_dists_ave$Interaction_status <- factor(pec_dists_ave$Interaction_status, 
                                           levels = c("interacting_all_species", 
                                                      "interacting_deca_only", 
                                                      "interacting_octbi_only", 
                                                      "not_interacting_Nany_species"))

bar_interactions_bins <- ggplot(pec_dists_ave, aes(x = Interaction_status, fill = Bin_Label_3)) +
  geom_bar(position = "dodge", stat = "count", alpha = 0.8, width = 0.7, color = "darkgrey") +
  labs(x = "Interaction status", y = "Count") +
  scale_fill_manual(values = bin_colors_3, 
                    labels = c(expression("< 5Mb apart in" ~ italic("P. maximus")),
                               expression("≥ 5Mb and < 15Mb apart in" ~ italic("P. maximus")),
                               expression("≥ 15Mb apart in" ~ italic("P. maximus")))) + theme_minimal() + 
  scale_x_discrete(labels = c("Interacting\nall species", "Interacting\ndecapodiformes","Interacting\nO. bimaculoides", "Interacting\nno species")) +
  theme(axis.title.x = element_text(size = 17.5),
        axis.title.y = element_text(size = 17.5),
        axis.text.x = element_text(size = 15.5, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15.5),
        legend.text = element_text(hjust = 0, size = 17),
        legend.title = element_blank()) + 
  geom_text(aes(label = ..count..), stat = "count", position = position_dodge(width = 0.7), size = 5.5, vjust = -0.5)


ggsave("bar_interactions_bins_10threshold_ave_with_sof.tiff", bar_interactions_bins, units = "in", dpi = 400, width = 16, height = 8)

