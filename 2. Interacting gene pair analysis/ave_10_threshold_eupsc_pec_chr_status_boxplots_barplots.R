# Classify different interacting gene pair categories, where gene pairs with 10 or more normalised read pairs as classified as interacting.
# Then plot boxplots of genomic distance between gene pairs in E. scolopes across interaction categories and colour by P. maximus chromosome status.
# Test for significant differences between distances for different P. maximus chromosome categories per interaction status using a pairwise Wilcoxon test with BH correction.
# Then plot barplots and stacked percentage barplots for the number of gene pairs each interaction category and P. maximus status.

# Clear workspace
rm(list = ls())

#nLoad libraries---
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyr)
library(tibble)

# Reading the output of synteny_by_topology_interactions.py script - interaction info based E. scolopes interactions 100 kb ----
eupsc_100k_obi_50k_pecten_chr <- read.delim("409493_intrachrom_allchrs_KR_100000_EUPvs50000_OBI_int_freq_PEC_chr_status_with_sof.txt")
head(eupsc_100k_obi_50k_pecten_chr)
dim(eupsc_100k_obi_50k_pecten_chr)

#Reading the output of eup_vs_obi_genom_dist_form.py  script -Based on E. scolopes interactions at 100 kb ----
eup_100k_obi_dist  <- read.delim("409493_intrachrom_allchrs_KR_100000_eupsc_octbi_genom_dist_sorted_merged_rm_dups.txt")

head(eup_100k_obi_dist)

#Merge datasets----
interactions_pecten_chr_genom_dist <- merge(eup_100k_obi_dist, eupsc_100k_obi_50k_pecten_chr, by=c(1))
head(interactions_pecten_chr_genom_dist)
dim(interactions_pecten_chr_genom_dist)

#Add column of interaction status based on your given threshold----
interactions_pecten_chr_genom_dist <- interactions_pecten_chr_genom_dist %>%
  mutate(
    Interaction_status = case_when(
      Interaction_frequency_spp1 >= 10 & Interaction_frequency_spp2 >= 10 & Interaction_frequency_Sepia >= 10 ~ 'interacting_all_species',
      Interaction_frequency_spp1 < 10 & Interaction_frequency_spp2 < 10 & Interaction_frequency_Sepia < 10 ~ 'not_interacting_any_species',
      Interaction_frequency_spp1 >= 10 & Interaction_frequency_spp2 < 10 & Interaction_frequency_Sepia >= 10 ~ 'interacting_deca_only',
      Interaction_frequency_spp1 < 10 & Interaction_frequency_spp2 >= 10 & Interaction_frequency_Sepia < 10 ~ 'interacting_octbi_only',
    )
  )

#Some checks----
# Count rows with defined conditions
rows_with_conditions <- interactions_pecten_chr_genom_dist %>%
  filter(!is.na(Interaction_status)) %>%
  nrow()

# Print the count
rows_with_conditions

# Count rows with defined conditions
rows_with_conditions <- sum(!is.na(interactions_pecten_chr_genom_dist$Interaction_status))

# Print the count
rows_with_conditions

# Remove rows without conditions (i.e. ones that interact in only Sepia or Euprymna)----
interactions_pecten_chr_genom_dist <- interactions_pecten_chr_genom_dist %>%
  filter(!is.na(Interaction_status))

# View the updated dataframe
head(interactions_pecten_chr_genom_dist)

# Average interaction frequencies for duplicate interaction frequencies in each category. This will unbias results caused by longer genes, but doing it for each category keeps cases of genes that may be on TAD borders etc.----
pec_dists_ave <- interactions_pecten_chr_genom_dist %>%
  group_by(Orth_pair_based_on_eupsc_interaction_matrix, Interaction_status, 
           Pecten_chromosome_status,
           Genomic_distance_eupsc_bp, Genomic_distance_octbi_bp) %>%
  summarize(
    Average_interaction_frequency_spp1 = mean(Interaction_frequency_spp1),
    Average_interaction_frequency_spp2 = mean(Interaction_frequency_spp2),
    Average_interaction_frequency_Sepia = mean(Interaction_frequency_Sepia)
  ) %>%
  ungroup()

head(pec_dists_ave)

# Save merged file for future analyses----
write.table(pec_dists_ave, "409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof.txt", append = FALSE, sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE, quote = FALSE)

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
  group_by(Pecten_chromosome_status) %>%
  summarize(count = n())
print(gene_counts)


# Subset to make different plots, because they don't all fit in the same one
pec_dists_ave_not_int <- subset(pec_dists_ave, subset=(Interaction_status=="not_interacting_any_species" | Interaction_status=="interacting_octbi_only" ))
pec_dists_ave_int <- subset(pec_dists_ave, subset=(Interaction_status=="interacting_all_species" | Interaction_status=="interacting_deca_only")) 

# Reordering the levels of Interaction_status so same pec chromosome is first in the plot. You have to also order Interaction_status too or it messes up the formatting.
pec_dists_ave_not_int <- transform(pec_dists_ave_not_int,
                                   Pecten_chromosome_status = factor(Pecten_chromosome_status, levels = c("same_pec_chrs", "diff_pec_chrs")),
                                   Interaction_status = factor(Interaction_status, levels = c("not_interacting_any_species", "interacting_octbi_only")))

pec_dists_ave_int <- transform(pec_dists_ave_int,
                               Pecten_chromosome_status = factor(Pecten_chromosome_status, levels = c("same_pec_chrs", "diff_pec_chrs")),
                               Interaction_status = factor(Interaction_status, levels = c("interacting_all_species", "interacting_deca_only")))

# Plotting the grouped boxplot
grouped_boxplot_not_int <- ggplot(pec_dists_ave_not_int, aes(x = Interaction_status, y = Genomic_distance_eupsc_bp, fill = Pecten_chromosome_status)) +
  geom_boxplot(outlier.size = 4, size = 2) +
  labs(x = NULL, y = expression(paste("Genomic distance in ", italic("E. scolopes"), " (bp)"))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), # Remove x-axis labels
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 50), # Adjust font size for y-axis labels
        axis.title.y = element_text(size = 50), # Adjust font size for y-axis title
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_fill_manual(values = c("#d62728", "#1f77b4"), # Adjust colors as needed
                    labels = c(expression(paste("Same ", italic("P. maximus"), " chromosome, n = 18997")), 
                               expression(paste("Different ", italic("P. maximus"), " chromosome, n = 40764")))) + # Use actual labels for statuses
  guides(fill = guide_legend(title = NULL)) + # Remove the legend title
  coord_cartesian(y = c(0,200000000))

# Save the plot
ggsave("grouped_box_>=10_pec_chr_status_ave_eupsc_dist_not_int_with_sof.tiff", 
       grouped_boxplot_not_int, 
       width = 15, 
       height = 30, 
       units = "in", 
       dpi = 300)


# Plotting the grouped boxplot
grouped_boxplot_int <- ggplot(pec_dists_ave_int, aes(x = Interaction_status, y = Genomic_distance_eupsc_bp, fill = Pecten_chromosome_status)) +
  geom_boxplot(outlier.size = 4, size = 2) +
  labs(x = NULL, y = expression(paste("Genomic distance in ", italic("E. scolopes"), " (bp)"))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), # Remove x-axis labels
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 50), # Adjust font size for y-axis labels
        axis.title.y = element_text(size = 50), # Adjust font size for y-axis title
        legend.position = "top", 
        legend.box = "horizontal") +
  scale_fill_manual(values = c("#d62728", "#1f77b4"), # Adjust colors as needed
                    labels = c(expression(paste("Same ", italic("P. maximus"), " chromosome, n = 18997")), 
                               expression(paste("Different ", italic("P. maximus"), " chromosome, n = 40764")))) + # Use actual labels for statuses
  guides(fill = guide_legend(title = NULL)) + # Remove the legend title
  coord_cartesian(y = c(0,7500000))

# Save the plot
ggsave("grouped_box_>=10_pec_chr_status_ave_eupsc_dist_int_with_sof.tiff", 
       grouped_boxplot_int, 
       width = 15, 
       height = 30, 
       units = "in", 
       dpi = 300)

# Perform pairwise Wilcoxon tests and BH correction
pairwise_results<- lapply(unique(pec_dists_ave$Interaction_status), function(status) {
  subset_data <- pec_dists_ave %>% filter(Interaction_status == status)
  
  # Perform the pairwise Wilcoxon test
  test_result <- pairwise.wilcox.test(subset_data$Genomic_distance_eupsc_bp, subset_data$Pecten_chromosome_status, p.adjust.method = "BH")
  
  # Extract p-values and format results
  result_df <- as.data.frame(test_result$p.value)
  
  # Convert the matrix to a long format data frame
  result_long <- result_df %>%
    tibble::rownames_to_column(var = "Group1") %>%
    pivot_longer(cols = -Group1, names_to = "Group2", values_to = "P_Value") %>%
    mutate(Interaction_status = status)
  
  result_long
})

# Combine results into a single data frame
pairwise_results_df <- bind_rows(pairwise_results)

print(pairwise_results_df)

write.table(pairwise_results_df, "pec_status_eupsc_dist_wilcox.txt", append = FALSE, sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Calculate medians for each group
median_results <- pec_dists_ave%>%
  group_by(Interaction_status, Pecten_chromosome_status) %>%
  summarize(Median_Genomic_Distance = median(Genomic_distance_eupsc_bp), .groups = 'drop')

# Print the medians
print(median_results)


# Grouped barplot----
pec_dists_ave$Pecten_chromosome_status <- factor(pec_dists_ave$Pecten_chromosome_status, 
                                                 levels = c("same_pec_chrs", "diff_pec_chrs"))  # Red first blue second
pec_dists_ave$Interaction_status <- factor(pec_dists_ave$Interaction_status, 
                                           levels = c("interacting_all_species", 
                                                      "interacting_deca_only", 
                                                      "interacting_octbi_only", 
                                                      "not_interacting_Nany_species"))

bar_interactions_chr_status_all <- ggplot(pec_dists_ave, aes(x = Interaction_status, fill = Pecten_chromosome_status)) +
  geom_bar(position = "dodge", stat = "count", alpha = 0.8, width = 0.7, color = "darkgrey") +
  labs(x = "Interaction status", y = "Count") +
  scale_fill_manual(values = c("same_pec_chrs" = "#E41A1C", "diff_pec_chrs" = "#377EB8")) +
  theme_minimal() + 
  scale_x_discrete(labels = c("Interacting\nall species", "Interacting\ndecapodiformes","Interacting\nO. bimaculoides", "Interacting\nno species")) +
  theme(axis.title.x = element_text(size=17.5), axis.title.y = element_text(size=17.5),
        axis.text.x = element_text(size = 15.5, angle = 45, hjust = 1), axis.text.y = element_text(size = 15.5),
        legend.text = element_text(hjust = 0, size=17),  # Align legend labels to the left
        legend.title = element_blank()) + 
  geom_text(aes(label = ..count..), stat = "count", position = position_dodge(width = 0.7), size = 5.5, vjust = -0.5)

ggsave("bar_interactions_chr_status_10threshold_ave_with_sof.tiff", bar_interactions_chr_status_all, units = "in", dpi = 300,  height = 8, width = 14)


# Stacked percentage barplot
# Calculate percentages
pec_dists_ave_percent <- pec_dists_ave %>%
  group_by(Interaction_status, Pecten_chromosome_status) %>%
  summarise(count = n()) %>%
  group_by(Interaction_status) %>%
  mutate(percent = count / sum(count) * 100)

# Order
pec_dists_ave_percent$Pecten_chromosome_status <- factor(
  pec_dists_ave_percent$Pecten_chromosome_status,
  levels = c("diff_pec_chrs","same_pec_chrs")  # Bottom to top in the stack
)

# Make the stacked barplot
bar_interactions_chr_status_stacked <- ggplot(pec_dists_ave_percent, aes(x = Interaction_status, y = percent, fill = Pecten_chromosome_status)) +
  geom_col(position = "stack", alpha = 0.9, color = "darkgrey", width = 0.7) +
  labs(x = "Interaction status", y = "Percentage") +
  scale_fill_manual(values = c("same_pec_chrs" = "#E41A1C", "diff_pec_chrs" = "#377EB8")) +
  scale_x_discrete(labels = c("Interacting\nacross coleoids", 
                              "Interacting\nDecapodiformes", 
                              "Interacting\nO. bimaculoides", 
                              "Interacting\nno coleoids")) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 17.5), 
        axis.title.y = element_text(size = 17.5),
        axis.text.x = element_text(size = 15.5, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 15.5),
        legend.text = element_text(hjust = 0, size = 17),  
        legend.title = element_blank())

print(bar_interactions_chr_status_stacked)

# Save the plot
ggsave("bar_interactions_chr_status_10threshold_ave_stacked_percent.tiff", 
       bar_interactions_chr_status_stacked, units = "in", dpi = 300, height = 8, width = 14)

