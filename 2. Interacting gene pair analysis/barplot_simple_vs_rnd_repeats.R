# Make bar plots of simple repeats vs. RND elements and perform Wilcoxon tests with BH correction to check if repeat type between interaction categories varies significantly.

# Clear workspace
rm(list = ls())

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

# Load data directly
all_summed_all <- read.delim("all_repeats_norm_eupsc100k/all_repeats_norm_eupsc100k_interacting_all_species.txt")
all_summed_none <- read.delim("all_repeats_norm_eupsc100k/all_repeats_norm_eupsc100k_not_interacting_any_species.txt")
all_summed_deca <- read.delim("all_repeats_norm_eupsc100k/all_repeats_norm_eupsc100k_interacting_deca_only.txt")
all_summed_octbi <- read.delim("all_repeats_norm_eupsc100k/all_repeats_norm_eupsc100k_interacting_octbi_only.txt")

# Add Interaction_status columns
all_summed_all$Interaction_status <- "interacting_all_species"
all_summed_none$Interaction_status <- "not_interacting_any_species"
all_summed_deca$Interaction_status <- "interacting_deca_only"
all_summed_octbi$Interaction_status <- "interacting_octbi_only"

# Combine data for later analysis
combined_data <- bind_rows(all_summed_all, all_summed_none, all_summed_deca, all_summed_octbi)
head(combined_data)

#Remove repeat counts of 1 because these rows skew the analysis 
combined_data <- combined_data[combined_data$TotalRepeatCount > 0, ]

# Categorise RepeatType as "simple Repeats" or "rnd Elements"
combined_data <- combined_data %>%
  mutate(RepeatCategory = ifelse(grepl("^rnd", RepeatType), "RND elements", "Simple repeats"))

# Summarise data by Interaction_status and RepeatCategory
summary_data <- combined_data %>%
  group_by(Interaction_status, RepeatCategory) %>%
  summarise(TotalNormCount = sum(NormCount), .groups = 'drop')

# Define the custom colour palette for interaction statuses
interaction_cols <- c(
  "interacting_all_species" = adjustcolor("#FF7878", alpha.f = 0.9), 
  "interacting_deca_only" = adjustcolor("palegreen", alpha.f = 0.9), 
  "interacting_octbi_only" = adjustcolor("mediumpurple1", alpha.f = 0.9), 
  "not_interacting_any_species" = adjustcolor("orange", alpha.f = 0.9)
)

# Create custom labels for each combination of interaction status and repeat category
summary_data$CustomLabel <- with(summary_data, paste0(RepeatCategory, " - ", Interaction_status))

# View bar plot of normalised repeat count per interaction category, where all repeat categories have the same colour
bar_plot  <- ggplot(summary_data, aes(x = RepeatCategory, y = TotalNormCount, fill = Interaction_status)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Interaction_status) +
  scale_fill_manual(values = interaction_cols) +
  theme_minimal() +
  labs(x = "Repeat category", y = "Total normalised count", fill = "Interaction status",
       title = "Normalised count of simple repeats vs. RND elements by interaction status in O. bimaculoides") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
print(bar_plot)

summary(combined_data$NormCount)


# Square root the data. This spreads out small numbers smoothly, unlike log which doesn't work so well for values under close to 0.
combined_data <- combined_data %>%
  mutate(SqrtNormCount = sqrt(NormCount))

# Perform Wilcoxon test for each interaction status
summary(combined_data$sqrtNormCount)
wilcox_results <- combined_data %>%
  group_by(Interaction_status) %>%
  summarise(p_value = wilcox.test(SqrtNormCount[RepeatCategory == "RND elements"],
                                  SqrtNormCount[RepeatCategory == "Simple repeats"])$p.value)

# Apply Benjamini-Hochberg correction
wilcox_results <- wilcox_results %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH"))

# Print results
print(wilcox_results)

# Summarise total count per category and convert to a percentage
enrichment_data <- combined_data %>%
  group_by(Interaction_status, RepeatCategory) %>%
  summarise(TotalCount = sum(NormCount), .groups = 'drop') %>%
  group_by(Interaction_status) %>%
  mutate(Percentage = TotalCount / sum(TotalCount) * 100)  # Convert to percentage

# Rename interaction status labels (remove underscores, improve readability)
enrichment_data$Interaction_status <- factor(enrichment_data$Interaction_status,
                                             levels = c("interacting_all_species",
                                                        "interacting_deca_only",
                                                        "interacting_octbi_only",
                                                        "not_interacting_any_species"),
                                             labels = c("Interacting all species",
                                                        "Interacting decapodiformes only",
                                                        "Interacting O. bimaculoides only",
                                                        "Not interacting any species"))


# Percentage stacked bar plot
stacked_barplot <- ggplot(enrichment_data, aes(x = Interaction_status, y = Percentage, fill = RepeatCategory)) +
  geom_bar(stat = "identity", position = "stack") +  
  scale_fill_manual(values = c("Simple repeats" = "#B2DFD8", "RND elements" = "#003F4F")) +
  theme_minimal() +
  labs(x = "Interaction status", y = "Percentage of repeat type",
       fill = "Repeat type",
       title = "Percentage of RND elements vs. simple repeats in each interaction status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


print(stacked_barplot)

# Save the plot
ggsave(stacked_barplot, file = "repeats_eupsc_with_sof.tiff", height = 4.5, width = 6, dpi = 200)





