# Plot density plot of mean insulation scores between gene pairs for different interaction categories and test for significant differences across categories using T test with BH correction

# Clear workplace 
rm(list = ls())

#Unload all previous packages if needed, and load others.
library(ggplot2)
library(dplyr)

# Reading the output of get_ave_ins_score_between_gene_pairs.py script----
eupsc_ins_score <- read.delim("409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof_eupsc_window350kb_ins_score.txt")
dim(eupsc_ins_score)
head(eupsc_ins_score)

# Subset the dataset to keep only the required columns so you can next remove duplicate rows, otherwise you will have duplicates in the file (same int freq in focal species, but different in the other)
eupsc_ins_score_with_int_status <- eupsc_ins_score[, c("Orth_pair_based_on_eupsc_interaction_matrix","Interaction_status","Pecten_chromosome_status", "Average_insulation_score_eupsc",  "Genomic_distance_eupsc_bp", "Average_interaction_frequency_spp1")]
eupsc_ins_score <- eupsc_ins_score[, c("Orth_pair_based_on_eupsc_interaction_matrix","Pecten_chromosome_status", "Average_insulation_score_eupsc", "Genomic_distance_eupsc_bp", "Average_interaction_frequency_spp1")]

# Remove Nans and make unique
eupsc_ins_score_with_int_status <- eupsc_ins_score_with_int_status[is.finite(eupsc_ins_score_with_int_status$Average_insulation_score_eupsc) & is.finite(eupsc_ins_score_with_int_status$Average_interaction_frequency_spp1), ]
eupsc_ins_score_with_int_status <- unique(eupsc_ins_score_with_int_status)

# Create the density plot----

# Define the custom colors for filling
interaction_cols <- c(
  "interacting_all_species" = adjustcolor("#FF7878", alpha.f = 0.7), 
  "interacting_deca_only" = adjustcolor("palegreen", alpha.f = 0.7), 
  "interacting_octbi_only" = adjustcolor("mediumpurple1", alpha.f = 0.7), 
  "not_interacting_any_species" = adjustcolor("orange", alpha.f = 0.7)
)

# Create the density plot with custom colors and italicized species name in x-axis label
ins_score_density <- ggplot(eupsc_ins_score_with_int_status, aes(x = Average_insulation_score_eupsc, fill = Interaction_status)) +
  geom_density(alpha = 0.7) +
  labs(title = "Density plot of average insulation score for each interaction status",
       x = expression(paste("Average insulation score in ", italic("E. scolopes"))),
       y = "Density",
       fill = "Interaction status") +
  scale_fill_manual(name = "Interaction status",
                    values = c(interaction_cols)) + 
  theme_minimal() + xlim(-2,2) +
  theme(
    plot.title = element_text(size = 20),      # Title text size
    axis.title.x = element_text(size = 20),    # X-axis title text size
    axis.title.y = element_text(size = 20),    # Y-axis title text size
    axis.text.x = element_text(size = 20),     # X-axis label text size
    axis.text.y = element_text(size = 20)      # Y-axis label text size
  )

# Save the plot
ggsave(ins_score_density, file = "ins_score_density_all_cats_eupsc_with_sof.tiff", height = 6, width = 10)


# Signficance tests on density distribution----
#  Perform pairwise t-tests (because it's a normal distribution) and compile results
ttest_results_eupsc <- combn(unique(eupsc_ins_score_with_int_status$Interaction_status), 2, FUN = function(pair) {
  group1 <- pair[1]
  group2 <- pair[2]
  
  group1_data <- eupsc_ins_score_with_int_status %>%
    filter(Interaction_status == group1) %>%
    pull(Average_insulation_score_eupsc)
  
  group2_data <- eupsc_ins_score_with_int_status %>%
    filter(Interaction_status == group2) %>%
    pull(Average_insulation_score_eupsc)
  
  t_test_result <- t.test(group1_data, group2_data)
  
  data.frame(
    Group1 = group1,
    Group2 = group2,
    Mean_Group1 = mean(group1_data),
    Mean_Group2 = mean(group2_data),
    Median_Group1 = median(group1_data),
    Median_Group2 = median(group2_data),
    P_Value = t_test_result$p.value
  )
}, simplify = FALSE)

# Combine results into a single data frame
ttest_results_df_eupsc <- do.call(rbind, ttest_results_eupsc)

# Apply Benjamini-Hochberg (BH) correction
ttest_results_df_eupsc$Adj_P_Value <- p.adjust(ttest_results_df_eupsc$P_Value, method = "BH")

# Print the results in a table and save
print(ttest_results_df_eupsc)

write.table(ttest_results_df_eupsc, "ns_core_eupsc_wilcox_with_sof.txt", append = FALSE, sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE, quote = FALSE)


