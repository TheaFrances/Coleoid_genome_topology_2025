# Make barplots of TAD boundary status in each interaction category using average insulation scores between gene pairs

# Clear workspace
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
both_spp_ins_score <- both_spp_ins_score[, c("Orth_pair_based_on_eupsc_interaction_matrix","Pecten_chromosome_status", "Average_insulation_score_eupsc", "Genomic_distance_eupsc_bp", "Average_insulation_score_octbi", "Genomic_distance_octbi_bp", "Average_insulation_score_sepof")]

# Remove Nans and make unique
both_spp_ins_score_with_int_status <- both_spp_ins_score_with_int_status[is.finite(both_spp_ins_score_with_int_status$Average_insulation_score_eupsc) & is.finite(both_spp_ins_score_with_int_status$Average_insulation_score_octbi) & is.finite(both_spp_ins_score$Average_insulation_score_sepof), ]
both_spp_ins_score_with_int_status <- unique(both_spp_ins_score_with_int_status)

both_spp_ins_score <- both_spp_ins_score[is.finite(both_spp_ins_score$Average_insulation_score_eupsc) & is.finite(both_spp_ins_score$Average_insulation_score_octbi) & is.finite(both_spp_ins_score$Average_insulation_score_sepof), ]
both_spp_ins_score <- unique(both_spp_ins_score)

# Create an empty data frame to store the counts
interaction_summary <- data.frame(
  Interaction_status = character(),
  Category = character(),
  Count = integer()
)

# Iterate over the interaction statuses and calculate counts
for (status in unique(both_spp_ins_score_with_int_status$Interaction_status)) {
  # Filter data for each interaction status
  subset_data <- both_spp_ins_score_with_int_status[both_spp_ins_score_with_int_status$Interaction_status == status, ]
  
  # Count for each category
  count_all_above <- sum(subset_data$Average_insulation_score_eupsc >= 0.2 & subset_data$Average_insulation_score_sepof >= 0.2 & subset_data$Average_insulation_score_octbi >= 0.2)
  count_both_below <- sum(subset_data$Average_insulation_score_eupsc <= -0.2 & subset_data$Average_insulation_score_sepof <= -0.2 & subset_data$Average_insulation_score_octbi <= -0.2)
  count_deca_above_octbi_below <- sum(subset_data$Average_insulation_score_eupsc >= 0.2 & subset_data$Average_insulation_score_sepof >= 0.2 & subset_data$Average_insulation_score_octbi <= -0.2)
  count_deca_below_octbi_above <- sum(subset_data$Average_insulation_score_eupsc <= -0.2 & subset_data$Average_insulation_score_sepof <= -0.2 & subset_data$Average_insulation_score_octbi >= 0.2)
  
  # Append the results to the interaction_summary dataframe
  interaction_summary <- rbind(interaction_summary,
                               data.frame(Interaction_status = status, Category = "All species >= 0.2", Count = count_all_above),
                               data.frame(Interaction_status = status, Category = "Deca >= 0.2 and Octbi =< -0.2", Count = count_deca_above_octbi_below),
                               data.frame(Interaction_status = status, Category = "Deca =< -0.2 and Octbi >= 0.2", Count = count_deca_below_octbi_above),
                               data.frame(Interaction_status = status, Category = "All species =< -0.2", Count = count_both_below)
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

# Grouped barplot
ggplot(interaction_summary, aes(x = Interaction_status, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Interaction status", y = "Count", fill = "Insulation score category") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12)
  )

# Define the custom colors for filling
interaction_cols <- c(
  "interacting_all_species" = adjustcolor("#FF7878", alpha.f = 0.9), 
  "interacting_deca_only" = adjustcolor("palegreen", alpha.f = 0.9), 
  "interacting_octbi_only" = adjustcolor("mediumpurple1", alpha.f = 0.9), 
  "not_interacting_any_species" = adjustcolor("orange", alpha.f = 0.9)
)

# Stacked barplot
int_status_by_tad_bar <- interaction_summary %>%
  group_by(Category) %>%
  mutate(Prop = Count / sum(Count)) %>%
  ggplot(aes(x = Category, y = Prop, fill = Interaction_status)) +  # Swap x and fill
  geom_bar(stat = "identity", position = "fill", width = 0.7) +  # Adjust bar width
  scale_y_continuous(labels = scales::percent_format()) +  # Convert y-axis to percentage
  scale_fill_manual(values = interaction_cols) +  # Use custom interaction colors
  labs(x = "Gene pair insulation score category", y = "Proportion", fill = "Interaction status") +  # Swap labels
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12)
  ) + 
  scale_x_discrete(labels = c("Within continous interaction\ndomain in all species", "Spanning insulating boundary\nin O. bimaculoides only","Spanning insulating boundary\nin decapodiformes only", "Spanning insualting\nboundary in all species"))

print(int_status_by_tad_bar)

# Save the plot
ggsave(int_status_by_tad_bar, file = "int_status_by_tad_bar_with_sof.tiff", height = 6, width = 8)

# Check the counts
interaction_summary_check <- both_spp_ins_score_with_int_status %>%
  group_by(Interaction_status) %>%
  summarise(
    count_all_above = sum(Average_insulation_score_eupsc >= 0.2 & 
                            Average_insulation_score_sepof >= 0.2 & 
                            Average_insulation_score_octbi >= 0.2),
    
    count_both_below = sum(Average_insulation_score_eupsc <= -0.2 & 
                             Average_insulation_score_sepof <= -0.2 & 
                             Average_insulation_score_octbi <= -0.2),
    
    count_deca_above_octbi_below = sum(Average_insulation_score_eupsc >= 0.2 & 
                                         Average_insulation_score_sepof >= 0.2 & 
                                         Average_insulation_score_octbi <= -0.2),
    
    count_deca_below_octbi_above = sum(Average_insulation_score_eupsc <= -0.2 & 
                                         Average_insulation_score_sepof <= -0.2 & 
                                         Average_insulation_score_octbi >= 0.2)
  )

print(interaction_summary_check)
