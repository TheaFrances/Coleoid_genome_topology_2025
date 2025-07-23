# Clear workspace
rm(list = ls())

#Load libraries----
library(ggplot2)
library(viridis)
library(MASS)
library(dplyr)
library(ggpubr)
library(scales)
library(data.table)

# Reading the output of add_sof_interactions.py script - interaction info based E. scolopes interactions 100 kb. Make some subsets of the data for use later. ----
interaction_freq_file <- read.delim("409493_intrachrom_allchrs_KR_100000_EUPvs50000_OBI_int_freq_PEC_chr_status_with_sof.txt")
head(interaction_freq_file)
dim(interaction_freq_file)


# Remove the octbi interaction frequency and then remove duplicate rows, otherwise df is too large, and you don't need to assess all 3 spp here.
interaction_freq_file <- interaction_freq_file %>%
  dplyr::select(-Interaction_frequency_spp2) %>% # Remove the column
  distinct()                              # Remove duplicate rows

# Log interaction frequency columns
interaction_freq_file <- interaction_freq_file %>%
  mutate(across(starts_with("Interaction_frequency"), ~ log(as.numeric(as.character(.)) + 0.01), .names = "log_{.col}"))

# Reading the output of eup_vs_obi_genom_dist_form.py  script - formatted version-----
# Note do not use filtered version with 'interaction categories' because the categories will not necessarily match interaction frequency when you merge (there is often more than one status per gene pair).
all_dists  <- read.delim("/Users/users/Desktop/Micro-C/topology_strength_analysis/eupsc_100k/409493_intrachrom_allchrs_KR_100000_interaction_freq_file_no_dups.txt")
 head(all_dists)
 head(interaction_freq_file)
 
# Merge datasets----
# First clean data to speed up process
all_dists[, 1] <- trimws(all_dists[, 1])
interaction_freq_file[, 1] <- trimws(interaction_freq_file[, 1])

# Convert data frames to data.tables and merge this way to keep values 
setDT(all_dists)
setDT(interaction_freq_file)

# Perform the join with allow.cartesian = TRUE
interactions_distance_merged <- merge(
  all_dists, 
  interaction_freq_file, 
  by.x = "Orth_pair_based_on_eupsc_interaction_matrix", 
  by.y = "Orth_interacting_gene_pair", 
  allow.cartesian = TRUE
)
head(interactions_distance_merged)

# Set distances to numeric and check maximum----
interactions_distance_merged$Genomic_distance_eupsc_bp <- as.numeric(interactions_distance_merged$Genomic_distance_eupsc_bp)
interactions_distance_merged$Genomic_distance_octbi_bp <- as.numeric(interactions_distance_merged$Genomic_distance_octbi_bp)
interactions_distance_merged$Genomic_distance_sepof_bp <- as.numeric(interactions_distance_merged$Genomic_distance_sepof_bp)

max(interactions_distance_merged$Genomic_distance_octbi_bp)
max(interactions_distance_merged$Genomic_distance_eupsc_bp)
max(interactions_distance_merged$Genomic_distance_sepof_bp)
max(interactions_distance_merged$Interaction_frequency_spp1)
max(interactions_distance_merged$Interaction_frequency_Sepia)


# Plots colouring different distances----
# Adding a new column 'dist_cat' (distance category) based on genomic distance between gene pairs
interactions_distance_merged_categories <- interactions_distance_merged%>%
  mutate(
    dist_cat = case_when(
      Genomic_distance_eupsc_bp >= 20000000 & Genomic_distance_sepof_bp >= 20000000 ~ 'over_20mb',
      Genomic_distance_eupsc_bp < 1000000 & Genomic_distance_sepof_bp < 10000000 ~ 'under_1mb',
      Genomic_distance_eupsc_bp >= 1000000 & Genomic_distance_eupsc_bp < 10000000 & Genomic_distance_sepof_bp >= 1000000 & Genomic_distance_sepof_bp < 10000000 ~ '1-10mb',
      Genomic_distance_eupsc_bp >= 10000000 & Genomic_distance_eupsc_bp < 20000000 & Genomic_distance_sepof_bp >= 10000000 & Genomic_distance_sepof_bp < 20000000 ~ '10-20mb',
    )
  )

# Examine plot of all genomic distances----
ggplot(interactions_distance_merged_categories, aes(x=log_Interaction_frequency_spp1, y=log_Interaction_frequency_Sepia,  color=dist_cat)) + 
  geom_point(size=1) + 
  scale_color_brewer(palette = "Set1", direction=-1)+
  theme(axis.title.x = element_text(size=12.5), axis.title.y = element_text(size=12.5),
        legend.text = element_text(hjust = 0),  # Align legend labels to the left
        legend.title = element_blank()) +  # Remove legend title
  ylab(expression(paste("Interaction frequency in ", italic("S. officinalis"), " (100 kb resolution)"))) +
  xlab(expression(paste(italic("E. scolopes "),  "interaction frequency (100 kb resolution)")))

# Plot without under 1Mb distances different colours and save----
interactions_distance_merged_over_1mb_categories  <- subset(interactions_distance_merged_categories, subset=!(dist_cat=="under_1mb"))

threshold_1_20mb <- ggplot(interactions_distance_merged_over_1mb_categories, aes(x=log_Interaction_frequency_spp1, y=log_Interaction_frequency_Sepia, color=dist_cat)) + 
  geom_point(size=1, alpha = 0.8) + 
  scale_color_manual(values = c("darkslategray3", "darkslategray4",  "darkslategray"),  labels = c("1-10Mb", "≥10-20Mb", "≥20Mb")) +  # Set custom color palette and key labels
  theme(axis.title.x = element_text(size=12.5), axis.title.y = element_text(size=12.5),
        legend.text = element_text(hjust = 0),
        legend.title = element_blank()) +
  ylab(expression(paste("Interaction frequency in ", italic("S. officinalis"), " (100  kb resolution, logged)"))) +
  xlab(expression(paste(italic("E. scolopes"), " interaction frequency (100 kb resolution, logged)"))) +
  xlim(-1.5,5)+ ylim(-1.5,5)

ggsave("interaction_threshold>=1Mb>=20Mb_eupsc_sepof_logged.tiff", threshold_1_20mb, units = "in", dpi = 200, width = 6.5, height = 5.6)

# Plot without under 10Mb distances and save----
interactions_distance_merged_over_10mb_categories <- subset(interactions_distance_merged_categories, dist_cat %in% c("over_20mb", "10-20mb"))

threshold_10_20mb <- ggplot(interactions_distance_merged_over_10mb_categories, aes(x=log_Interaction_frequency_spp1, y=log_Interaction_frequency_Sepia, color=dist_cat)) + 
  geom_point(size=1, alpha = 0.8) + 
  scale_color_manual(values = c("darkslategray4", "darkslategray"), labels = c("≥10-20Mb", "≥20Mb")) +  # Set custom color palette and key labels
  theme(axis.title.x = element_text(size=12.5), axis.title.y = element_text(size=12.5),
        legend.text = element_text(hjust = 0),
        legend.title = element_blank()) +
  ylab(expression(paste("Interaction frequency in ", italic("S. officinalis"), " (100 kb resolution, logged)"))) +
  xlab(expression(paste(italic("E. scolopes"), " interaction frequency (100 kb resolution, logged)"))) + xlim(-1.5,5)+ ylim(-1.5,5)

ggsave("interaction_threshold>=10Mb_eupsc_sepof_logged.tiff", threshold_10_20mb, units = "in", dpi = 200, width = 6.5, height = 5.6)


