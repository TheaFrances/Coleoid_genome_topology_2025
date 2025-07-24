#Make proportional barplots of TAD boundary status in each P. maximus chromosome status or intrachromosomal P. maximus distance category, using average insulation scores between gene pairs

#Clear workspace
rm(list = ls())

library(dplyr)
library(ggplot2)
library(scales)  # For percentage formatting

# Reading the output of get_ave_ins_score_between_gene_pairs.py script for all species----
eupsc_ins_score <- read.delim("409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof_eupsc_window350kb_ins_score.txt")
octbi_ins_score <- read.delim("409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof_octbi_window350kb_ins_score.txt")
sepof_ins_score <- read.delim("409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof_sepof_window350kb_ins_score.txt")

# Merge files
two_spp_ins_score <- merge(eupsc_ins_score,octbi_ins_score)
both_spp_ins_score <- merge(two_spp_ins_score,sepof_ins_score)
dim(both_spp_ins_score)
head(both_spp_ins_score)

# Subset the dataset to keep only the required columns so you can next remove duplicate rows, otherwise you will have duplicates in the file (same int freq in focal species, but different in the other)
both_spp_ins_score_with_int_status <- both_spp_ins_score[, c("Orth_pair_based_on_eupsc_interaction_matrix","Interaction_status","Pecten_chromosome_status", "Average_insulation_score_eupsc",  "Genomic_distance_eupsc_bp", "Average_insulation_score_octbi",  "Genomic_distance_octbi_bp", "Average_insulation_score_sepof")]

# Remove Nans and make unique
both_spp_ins_score_with_int_status <- both_spp_ins_score_with_int_status[is.finite(both_spp_ins_score_with_int_status$Average_insulation_score_eupsc) & is.finite(both_spp_ins_score_with_int_status$Average_insulation_score_octbi) & is.finite(both_spp_ins_score$Average_insulation_score_sepof), ]
both_spp_ins_score_with_int_status <- unique(both_spp_ins_score_with_int_status)

# Create an empty data frame to store the counts
interaction_summary <- data.frame(
  Interaction_status = character(),
  Category = character(),
  Count = integer()
)

# Iterate over the interaction statuses and calculate counts
for (status in unique(both_spp_ins_score_with_int_status$Pecten_chromosome_status)) {
  # Filter data for each interaction status
  subset_data <- both_spp_ins_score_with_int_status[both_spp_ins_score_with_int_status$Pecten_chromosome_status == status, ]
  
  # Count for each category
  count_all_above <- sum(subset_data$Average_insulation_score_eupsc >= 0.2 & subset_data$Average_insulation_score_sepof >= 0.2 & subset_data$Average_insulation_score_octbi >= 0.2)
  count_both_below <- sum(subset_data$Average_insulation_score_eupsc <= -0.2 & subset_data$Average_insulation_score_sepof <= -0.2 & subset_data$Average_insulation_score_octbi <= -0.2)
  count_deca_above_octbi_below <- sum(subset_data$Average_insulation_score_eupsc >= 0.2 & subset_data$Average_insulation_score_sepof >= 0.2 & subset_data$Average_insulation_score_octbi <= -0.2)
  count_deca_below_octbi_above <- sum(subset_data$Average_insulation_score_eupsc <= -0.2 & subset_data$Average_insulation_score_sepof <= -0.2 & subset_data$Average_insulation_score_octbi >= 0.2)
  
  # Append the results to the interaction_summary dataframe
  interaction_summary <- rbind(interaction_summary,
                               data.frame(Pecten_chromosome_status = status, Category = "All species >= 0.2", Count = count_all_above),
                               data.frame(Pecten_chromosome_status = status, Category = "Deca >= 0.2 and Octbi =< -0.2", Count = count_deca_above_octbi_below),
                               data.frame(Pecten_chromosome_status = status, Category = "Deca =< -0.2 and Octbi >= 0.2", Count = count_deca_below_octbi_above),
                               data.frame(Pecten_chromosome_status = status, Category = "All species =< -0.2", Count = count_both_below)
  )
}

# View the interaction_summary
print(interaction_summary)

# Reorder the Category levels
interaction_summary$Category <- factor(interaction_summary$Category, 
                                       levels = c("All species >= 0.2", 
                                                  "Deca >= 0.2 and Octbi =< -0.2", 
                                                  "Deca =< -0.2 and Octbi >= 0.2", 
                                                  "All species =< -0.2"))



# Define the custom colors for filling
pec_chr_stat_cols <- c(
  "diff_pec_chrs" = adjustcolor("#377EB8", alpha.f = 0.9), 
  "same_pec_chrs" = adjustcolor("#E41A1C", alpha.f = 0.9))


interaction_summary$Pecten_chromosome_status <- factor(interaction_summary$Pecten_chromosome_status)

# Stacked barplot
pec_chr_status_by_tad_bar <- interaction_summary %>%
  group_by(Category) %>%
  mutate(Prop = Count / sum(Count)) %>%
  ggplot(aes(x = Category, y = Prop, fill = Pecten_chromosome_status)) +  # Swap x and fill
  geom_bar(stat = "identity", position = "fill", width = 0.7) +  # Adjust bar width
  scale_y_continuous(labels = scales::percent_format()) +  # Convert y-axis to percentage
  scale_fill_manual(values = pec_chr_stat_cols) +  # Use custom interaction colors
  labs(x = "Gene pair insulation score category", y = "Proportion", fill = "Pecten_chromosome_status") +  # Swap labels
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12)
  ) + 
  scale_x_discrete(labels = c("Within continous interaction\ndomain in all species", "Spanning insulating boundary\nin O. bimaculoides only","Spanning insulating boundary\nin decapodiformes only", "Spanning insualting\nboundary in all species"))

print(pec_chr_status_by_tad_bar)

# Save the plot
ggsave(pec_chr_status_by_tad_bar, file = "pec_chr_status_by_tad_bar_with_sof.tiff", height = 6, width = 8)


#Now, instead of P. maximus chromosome status, split into P. maxmimus distance categories ####----

# Clear workspace
rm(list = ls())

# Reading the output of get_ave_ins_score_between_gene_pairs.py script - E. scolopes interactions 100Kb based on eupsc 100kb interaction matrix----
eupsc_ins_score <- read.delim("409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof_eupsc_window350kb_ins_score.txt")
octbi_ins_score <- read.delim("409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof_octbi_window350kb_ins_score.txt")
sepof_ins_score <- read.delim("409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof_sepof_window350kb_ins_score.txt")

#Merge files
two_spp_ins_score <- merge(eupsc_ins_score,octbi_ins_score)
both_spp_ins_score <- merge(two_spp_ins_score,sepof_ins_score)
dim(both_spp_ins_score)
head(both_spp_ins_score)

#Upload P. max distances file
pec_bins <- read.delim("/Users/users/Desktop/Micro-C/topology_strength_analysis/eupsc_100k/409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_pec_bins_with_sof.txt")

#Keep first and last column only
pec_bins <- pec_bins %>% select(1, last_col())

#Merge
spp_ins_score_bins <- merge(both_spp_ins_score, pec_bins)
head(spp_ins_score_bins)

#Subset the dataset to keep only the required columns so you can next remove duplicate rows, otherwise you will have duplicates in the file (same int freq in focal species, but different in the other)
spp_ins_score_bins_with_int_status <- spp_ins_score_bins[, c("Orth_pair_based_on_eupsc_interaction_matrix","Interaction_status","Bin_Label_3", "Average_insulation_score_eupsc",  "Genomic_distance_eupsc_bp", "Average_insulation_score_octbi",  "Genomic_distance_octbi_bp", "Average_insulation_score_sepof")]

#Remove Nans and make unique
spp_ins_score_bins_with_int_status <- spp_ins_score_bins_with_int_status[is.finite(spp_ins_score_bins_with_int_status$Average_insulation_score_eupsc) & is.finite(spp_ins_score_bins_with_int_status$Average_insulation_score_octbi) & is.finite(spp_ins_score_bins$Average_insulation_score_sepof), ]
spp_ins_score_bins_with_int_status <- unique(spp_ins_score_bins_with_int_status)


# Create an empty data frame to store the counts
interaction_summary <- data.frame(
  Interaction_status = character(),
  Category = character(),
  Count = integer()
)

# Iterate over the interaction statuses and calculate counts
for (status in unique(spp_ins_score_bins_with_int_status$Bin_Label_3)) {
  # Filter data for each interaction status
  subset_data <- spp_ins_score_bins_with_int_status[spp_ins_score_bins_with_int_status$Bin_Label_3 == status, ]
  
  # Count for each category
  count_all_above <- sum(subset_data$Average_insulation_score_eupsc >= 0.2 & subset_data$Average_insulation_score_sepof >= 0.2 & subset_data$Average_insulation_score_octbi >= 0.2)
  count_both_below <- sum(subset_data$Average_insulation_score_eupsc <= -0.2 & subset_data$Average_insulation_score_sepof <= -0.2 & subset_data$Average_insulation_score_octbi <= -0.2)
  count_deca_above_octbi_below <- sum(subset_data$Average_insulation_score_eupsc >= 0.2 & subset_data$Average_insulation_score_sepof >= 0.2 & subset_data$Average_insulation_score_octbi <= -0.2)
  count_deca_below_octbi_above <- sum(subset_data$Average_insulation_score_eupsc <= -0.2 & subset_data$Average_insulation_score_sepof <= -0.2 & subset_data$Average_insulation_score_octbi >= 0.2)
  
  # Append the results to the interaction_summary dataframe
  interaction_summary <- rbind(interaction_summary,
                               data.frame(Bin_Label_3 = status, Category = "All species >= 0.2", Count = count_all_above),
                               data.frame(Bin_Label_3 = status, Category = "Deca >= 0.2 and Octbi =< -0.2", Count = count_deca_above_octbi_below),
                               data.frame(Bin_Label_3 = status, Category = "Deca =< -0.2 and Octbi >= 0.2", Count = count_deca_below_octbi_above),
                               data.frame(Bin_Label_3 = status, Category = "All species =< -0.2", Count = count_both_below)
  )
}

# View the interaction_summary
print(interaction_summary)

# Reorder the category levels
interaction_summary$Category <- factor(interaction_summary$Category, 
                                       levels = c("All species >= 0.2", 
                                                  "Deca >= 0.2 and Octbi =< -0.2", 
                                                  "Deca =< -0.2 and Octbi >= 0.2", 
                                                  "All species =< -0.2"))

# Define the custom colors for filling
pec_bin_cols <- c(
  "bin1" = adjustcolor("yellow3", alpha.f = 0.9), 
  "bin2" = adjustcolor("cyan4", alpha.f = 0.9),
  "bin3" = adjustcolor("darkslateblue", alpha.f = 0.9))

interaction_summary$Bin_Label_3 <- factor(interaction_summary$Bin_Label_3)

# Stacked barplot
int_status_by_tad_bar_pec_bin_status <- interaction_summary %>%
  group_by(Category) %>%
  mutate(Prop = Count / sum(Count)) %>%
  ggplot(aes(x = Category, y = Prop, fill = Bin_Label_3)) +  # Swap x and fill
  geom_bar(stat = "identity", position = "fill", width = 0.8) +  # Adjust bar width
  scale_y_continuous(labels = scales::percent_format()) +  # Convert y-axis to percentage
  scale_fill_manual(values = pec_bin_cols) +  # Use custom interaction colors
  labs(x = "Gene pair insulation score category", y = "Proportion", fill = "Bin_Label_3") +  # Swap labels
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 15)
  ) + 
  scale_x_discrete(labels = c("Within continous interaction\ndomain in all species", "Spanning insulating boundary\nin O. bimaculoides only","Spanning insulating boundary\nin decapodiformes only", "Spanning insualting\nboundary in all species"))

print(int_status_by_tad_bar_pec_bin_status)

ggsave(int_status_by_tad_bar_pec_bin_status, file = "pec_bin_status_by_tad_bar_with_sof.tiff", height = 7, width = 8)


