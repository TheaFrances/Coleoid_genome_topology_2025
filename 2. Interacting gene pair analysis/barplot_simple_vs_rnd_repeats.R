# Make bar plots of simple repeats vs. RND elements

# Clear workspace
rm(list = ls())

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

# Load data directly
all_summed_all  <- read.delim("all_repeats_norm_eupsc100k/all_repeats_norm_eupsc100k_interacting_all_species.txt")
all_summed_none <- read.delim("all_repeats_norm_eupsc100k/all_repeats_norm_eupsc100k_not_interacting_any_species.txt")
all_summed_deca <- read.delim("all_repeats_norm_eupsc100k/all_repeats_norm_eupsc100k_interacting_deca_only.txt")
all_summed_octbi<- read.delim("all_repeats_norm_eupsc100k/all_repeats_norm_eupsc100k_interacting_octbi_only.txt")

# Add Interaction_status columns
all_summed_all$Interaction_status   <- "interacting_all_species"
all_summed_none$Interaction_status  <- "not_interacting_any_species"
all_summed_deca$Interaction_status  <- "interacting_deca_only"
all_summed_octbi$Interaction_status <- "interacting_octbi_only"

# Combine data
combined_data <- bind_rows(all_summed_all, all_summed_none, all_summed_deca, all_summed_octbi)

# Remove rows with zero total repeat count (avoids empty bins)
combined_data <- combined_data[combined_data$TotalRepeatCount > 0, ]

# Categorise RepeatType as "Simple repeats" or "RND elements"
# (Everything not starting with "rnd" will be treated as "Simple repeats" here)
combined_data <- combined_data %>%
  mutate(RepeatCategory = ifelse(grepl("^rnd", RepeatType, ignore.case = TRUE),
                                 "RND elements", "Simple repeats"))

# --- Panel 1: Total normalised counts (faceted) ---
summary_data <- combined_data %>%
  group_by(Interaction_status, RepeatCategory) %>%
  summarise(TotalNormCount = sum(NormCount), .groups = 'drop')

# Custom colours for interaction statuses
interaction_cols <- c(
  "interacting_all_species"   = adjustcolor("#FF7878", alpha.f = 0.9),
  "interacting_deca_only"     = adjustcolor("palegreen", alpha.f = 0.9),
  "interacting_octbi_only"    = adjustcolor("mediumpurple1", alpha.f = 0.9),
  "not_interacting_any_species" = adjustcolor("orange", alpha.f = 0.9)
)

# Bar plot of normalised repeat counts
bar_plot <- ggplot(summary_data, aes(x = RepeatCategory, y = TotalNormCount, fill = Interaction_status)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Interaction_status) +
  scale_fill_manual(values = interaction_cols) +
  theme_minimal() +
  labs(x = "Repeat category", y = "Total normalised count", fill = "Interaction status",
       title = "Normalised count of simple repeats vs. RND elements by interaction status in O. bimaculoides") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(bar_plot)

# --- Panel 2: Percentages within each interaction status (stacked barplot) ---
enrichment_data <- combined_data %>%
  group_by(Interaction_status, RepeatCategory) %>%
  summarise(TotalCount = sum(NormCount), .groups = 'drop') %>%
  group_by(Interaction_status) %>%
  mutate(Percentage = TotalCount / sum(TotalCount) * 100)

# Rename interaction status labels for readability
enrichment_data$Interaction_status <- factor(
  enrichment_data$Interaction_status,
  levels = c("interacting_all_species",
             "interacting_deca_only",
             "interacting_octbi_only",
             "not_interacting_any_species"),
  labels = c("Interacting all species",
             "Interacting decapodiformes only",
             "Interacting O. bimaculoides only",
             "Not interacting any species")
)

# Stacked percentage bar plot
stacked_barplot <- ggplot(enrichment_data, aes(x = Interaction_status, y = Percentage, fill = RepeatCategory)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Simple repeats" = "#B2DFD8", "RND elements" = "#003F4F")) +
  theme_minimal() +
  labs(x = "Interaction status", y = "Percentage of repeat type",
       fill = "Repeat type",
       title = "Percentage of RND elements vs. simple repeats in each interaction status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(stacked_barplot)

# Save the stacked bar plot (fixing ggsave usage)
ggsave(filename = "repeats_eupsc_with_sof.tiff", plot = stacked_barplot, height = 4.5, width = 6, dpi = 200)
