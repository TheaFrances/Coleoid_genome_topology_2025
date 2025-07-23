rm(list = ls())

#lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE) #. Use if there's an error with ggplot.

#Load libraries part 1----
library(ggplot2)
library(viridis)
library(MASS)
library(dplyr)
library(ggpubr)
library(scales)


#Reading the output of add_sof_interactions.py script - interaction info based E. scolopes interactions 100 kb. Make some subsets of the data for use later. ----
eupsc_100k_obi_100k_pecten_chr <- read.delim("/Users/users/Desktop/Micro-C/topology_strength_analysis/eupsc_100k/409493_intrachrom_allchrs_KR_100000_EUPvs50000_OBI_int_freq_PEC_chr_status_with_sof.txt")
head(eupsc_100k_obi_100k_pecten_chr)
dim(eupsc_100k_obi_100k_pecten_chr)


# Remove the octbi interaction frequency and then remove duplicate rows, otherwise df is too large, and you don't need to assess all 3 spp here.
eupsc_100k_obi_100k_pecten_chr <- eupsc_100k_obi_100k_pecten_chr %>%
  dplyr::select(-Interaction_frequency_spp2) %>% # Remove the column
  distinct()                              # Remove duplicate rows

#Log interaction frequency columns
eupsc_100k_obi_100k_pecten_chr <- eupsc_100k_obi_100k_pecten_chr %>%
  mutate(across(starts_with("Interaction_frequency"), ~ log(as.numeric(as.character(.)) + 0.01), .names = "log_{.col}"))

#Reading the output of eup_vs_obi_genom_dist.py  script - formatted version- Based on eupsc interactions at 100----
#Note do not use filtered version with 'interaction categories' because the categories will not necessarily match interaction frequency when you merge (there is often more than one status per gene pair).
all_dists  <- read.delim("/Users/users/Desktop/Micro-C/topology_strength_analysis/eupsc_100k/409493_intrachrom_allchrs_KR_100000_eupsc_sepof_dists_no_dups.txt")
 head(all_dists)
 head(eupsc_100k_obi_100k_pecten_chr)
 
#Merge datasets----

#First clean data to speed up process
all_dists[, 1] <- trimws(all_dists[, 1])
eupsc_100k_obi_100k_pecten_chr[, 1] <- trimws(eupsc_100k_obi_100k_pecten_chr[, 1])

# Convert data frames to data.tables and merge this way to keep values 
library(data.table)
setDT(all_dists)
setDT(eupsc_100k_obi_100k_pecten_chr)

# Perform the join with allow.cartesian = TRUE
interactions_pecten_chr_genom_dist <- merge(
  all_dists, 
  eupsc_100k_obi_100k_pecten_chr, 
  by.x = "Orth_pair_based_on_eupsc_interaction_matrix", 
  by.y = "Orth_interacting_gene_pair", 
  allow.cartesian = TRUE
)
head(interactions_pecten_chr_genom_dist)

#Set distances to numeric and get max
interactions_pecten_chr_genom_dist$Genomic_distance_eupsc_bp <- as.numeric(interactions_pecten_chr_genom_dist$Genomic_distance_eupsc_bp)
interactions_pecten_chr_genom_dist$Genomic_distance_octbi_bp <- as.numeric(interactions_pecten_chr_genom_dist$Genomic_distance_octbi_bp)
interactions_pecten_chr_genom_dist$Genomic_distance_sepof_bp <- as.numeric(interactions_pecten_chr_genom_dist$Genomic_distance_sepof_bp)

max(interactions_pecten_chr_genom_dist$Genomic_distance_octbi_bp)
max(interactions_pecten_chr_genom_dist$Genomic_distance_eupsc_bp)
max(interactions_pecten_chr_genom_dist$Genomic_distance_sepof_bp)
max(interactions_pecten_chr_genom_dist$Interaction_frequency_spp1)
max(interactions_pecten_chr_genom_dist$Interaction_frequency_Sepia)


#Coloring different distances
# Adding a new column 'interaction_status' based on conditions
interactions_pecten_chr_genom_dist_categories <- interactions_pecten_chr_genom_dist %>%
  mutate(
    dist_cat = case_when(
      Genomic_distance_eupsc_bp >= 20000000 & Genomic_distance_sepof_bp >= 20000000 ~ 'over_20mb',
      Genomic_distance_eupsc_bp < 1000000 & Genomic_distance_sepof_bp < 10000000 ~ 'under_1mb',
      Genomic_distance_eupsc_bp >= 1000000 & Genomic_distance_eupsc_bp < 10000000 & Genomic_distance_sepof_bp >= 1000000 & Genomic_distance_sepof_bp < 10000000 ~ '1-10mb',
      Genomic_distance_eupsc_bp >= 10000000 & Genomic_distance_eupsc_bp < 20000000 & Genomic_distance_sepof_bp >= 10000000 & Genomic_distance_sepof_bp < 20000000 ~ '10-20mb',
    )
  )


ggplot(interactions_pecten_chr_genom_dist_categories, aes(x=log_Interaction_frequency_spp1, y=log_Interaction_frequency_Sepia,  color=dist_cat)) + 
  geom_point(size=1) + 
  scale_color_brewer(palette = "Set1", direction=-1)+
  theme(axis.title.x = element_text(size=12.5), axis.title.y = element_text(size=12.5),
        legend.text = element_text(hjust = 0),  # Align legend labels to the left
        legend.title = element_blank()) +  # Remove legend title
  ylab(expression(paste("Interaction frequency in ", italic("S. officinalis"), " (100 kb resolution)"))) +
  xlab(expression(paste(italic("E. scolopes "),  "interaction frequency (100 kb resolution)")))

#Without under 1Mb different colours
interactions_pecten_chr_genom_dist_over_1mb_categories  <- subset(interactions_pecten_chr_genom_dist_categories, subset=!(dist_cat=="under_1mb"))

threshold_1_20mb <- ggplot(interactions_pecten_chr_genom_dist_over_1mb_categories, aes(x=log_Interaction_frequency_spp1, y=log_Interaction_frequency_Sepia, color=dist_cat)) + 
  geom_point(size=1, alpha = 0.8) + 
  scale_color_manual(values = c("darkslategray3", "darkslategray4",  "darkslategray"),  labels = c("1-10Mb", "≥10-20Mb", "≥20Mb")) +  # Set custom color palette and key labels
  theme(axis.title.x = element_text(size=12.5), axis.title.y = element_text(size=12.5),
        legend.text = element_text(hjust = 0),
        legend.title = element_blank()) +
  ylab(expression(paste("Interaction frequency in ", italic("S. officinalis"), " (100  kb resolution, logged)"))) +
  xlab(expression(paste(italic("E. scolopes"), " interaction frequency (100 kb resolution, logged)"))) +
  xlim(-1.5,5)+ ylim(-1.5,5)

ggsave("/Users/users/Desktop/Micro-C/figs_for_paper/Figure3/interaction_threshold>=1Mb>=20Mb_with_sepof_logged.tiff", threshold_1_20mb, units = "in", dpi = 400, width = 6.5, height = 5.6)

#Without under 10Mb different colours
interactions_pecten_chr_genom_dist_over_10mb_categories <- subset(interactions_pecten_chr_genom_dist_categories, dist_cat %in% c("over_20mb", "10-20mb"))

# Plot scatterplot with custom colors and custom key labels
threshold_10_20mb <- ggplot(interactions_pecten_chr_genom_dist_over_10mb_categories, aes(x=log_Interaction_frequency_spp1, y=log_Interaction_frequency_Sepia, color=dist_cat)) + 
  geom_point(size=1, alpha = 0.8) + 
  scale_color_manual(values = c("darkslategray4", "darkslategray"), labels = c("≥10-20Mb", "≥20Mb")) +  # Set custom color palette and key labels
  theme(axis.title.x = element_text(size=12.5), axis.title.y = element_text(size=12.5),
        legend.text = element_text(hjust = 0),
        legend.title = element_blank()) +
  ylab(expression(paste("Interaction frequency in ", italic("S. officinalis"), " (100 kb resolution, logged)"))) +
  xlab(expression(paste(italic("E. scolopes"), " interaction frequency (100 kb resolution, logged)"))) + xlim(-1.5,5)+ ylim(-1.5,5)

ggsave("/Users/users/Desktop/Micro-C/figs_for_paper/Figure3/interaction_threshold>=10Mb_with_sepof_logged.tiff", threshold_10_20mb, units = "in", dpi = 400, width = 6.5, height = 5.6)


#Different way
#Divide interaction strength by distance and plot the ratios for esc (on x axis) and for obi (y axis), so ratios instead of interaction frequency/strength
#Here, small genomic distance/large interaction frequency will be outliers in the graph.

interactions_pecten_chr_genom_dist_ratios <- interactions_pecten_chr_genom_dist %>%
  mutate(Ratio_freq1_eupsc = Interaction_frequency_spp1 / Genomic_distance_eupsc_bp,
         Ratio_freq2_octbi = Interaction_frequency_Sepia / Genomic_distance_octbi_bp)

head(interactions_pecten_chr_genom_dist_ratios)

ggplot(interactions_pecten_chr_genom_dist_ratios, aes(x=Ratio_freq1_eupsc, y=Ratio_freq2_octbi)) + 
  geom_point(size=1) + 
  scale_color_brewer(palette = "Set1", direction=-1)+
  theme(axis.title.x = element_text(size=12.5), axis.title.y = element_text(size=12.5),
        legend.text = element_text(hjust = 0),  # Align legend labels to the left
        legend.title = element_blank()) +  # Remove legend title
  ylab(expression(paste("Interaction frequency in ", italic("S. officinalis"), " (100 kb resolution, logged)"))) +
  xlab(expression(paste(italic("E. scolopes "),  "interaction frequency (100 kb resolution, logged)")))

ggplot(interactions_pecten_chr_genom_dist_ratios, aes(x = Ratio_freq1_eupsc, y = Ratio_freq2_octbi)) +
  geom_point(size = 1) +
  labs(x = "Ratio of Interaction Frequency 1 (Spp1) to Genomic Distance (Eupsc)",
       y = "Ratio of Interaction Frequency 2 (Sepia) to Genomic Distance (Octbi)",
       title = "Scatterplot of Ratios") +
  theme_minimal() +
  coord_cartesian(xlim = c(0, 5), ylim = c(0, 5))

head(interactions_pecten_chr_genom_dist_ratios)
