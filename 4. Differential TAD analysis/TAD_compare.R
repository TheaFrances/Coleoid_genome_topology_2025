rm(list = ls())

# Load required packages
library(ggplot2)
library(TADCompare)
library(dplyr)
library(tidyr)
library(stringr)

# Get list of all matching input files
file_list <- Sys.glob("/Users/users/Desktop/Micro-C/scolopes/microc*/KR_norm/*/dumped_100k/*_KR_100000_cleaned.dumped.hic")

# Extract unique tissues from file paths
tissues <- unique(str_extract(file_list, "(?<=KR_norm/)[^/]+"))

# Generate all pairwise comparisons
comparisons <- expand.grid(Tissue1 = tissues, Tissue2 = tissues) %>%
  filter(Tissue1 != Tissue2) %>%  # Remove self-comparisons (e.g., "cb vs cb")
  rowwise() %>%
  mutate(Comparison = paste(sort(c(Tissue1, Tissue2)), collapse = " vs ")) %>%  # Standardize order
  distinct(Comparison, .keep_all = TRUE) %>%  # Remove duplicate reverse comparisons
  select(-Comparison)  # Remove extra column after filtering

# Initialize an empty data frame to store all results
all_plot_data <- data.frame()

# Assuming file_names contains full paths, extract just the chromosome names
chromosomes <- gsub(".*_(Lachesis_group[0-9]+__[0-9]+_contigs__length_[0-9]+).*", "\\1", file_list)

# Loop through each tissue comparison
for (i in 1:nrow(comparisons)) {
  tissue1 <- comparisons$Tissue1[i]
  tissue2 <- comparisons$Tissue2[i]
  print(paste("Processing comparison:", tissue1, "vs", tissue2))
  
  # Find file paths dynamically for both tissues
  file1 <- Sys.glob(paste0("/Users/users/Desktop/Micro-C/scolopes/microc*/KR_norm/", tissue1, "/dumped_100k/*_KR_100000_cleaned.dumped.hic"))
  file2 <- Sys.glob(paste0("/Users/users/Desktop/Micro-C/scolopes/microc*/KR_norm/", tissue2, "/dumped_100k/*_KR_100000_cleaned.dumped.hic"))
  
  # Skip if any tissue has no valid files
  if (length(file1) == 0 | length(file2) == 0) {
    next
  }
  
  # Loop through each chromosome
  for (chr in chromosomes) {
    chr_file1 <- grep(chr, file1, value = TRUE)
    chr_file2 <- grep(chr, file2, value = TRUE)
    
    # Skip missing chromosome files
    if (length(chr_file1) == 0 | length(chr_file2) == 0) {
      next
    }
    
    # Load the contact matrices
    contact_matrix1 <- read.table(chr_file1, header = FALSE)
    contact_matrix2 <- read.table(chr_file2, header = FALSE)
    
    # Find max bin position
    max_bin1 <- max(contact_matrix1[, 1:2])
    max_bin2 <- max(contact_matrix2[, 1:2])
    resolution <- 100000  # 100kb resolution
    window_size <- 10  # 10-bin window (1mb)
    n_bins <- ceiling(max(max_bin1, max_bin2) / resolution)
    
    # Create a vector of bin positions
    bin_positions <- seq(0, by = resolution, length.out = n_bins)
    
    # Create empty matrices
    dense_matrix1 <- matrix(0, nrow = n_bins, ncol = n_bins)
    rownames(dense_matrix1) <- bin_positions
    colnames(dense_matrix1) <- bin_positions
    
    dense_matrix2 <- matrix(0, nrow = n_bins, ncol = n_bins)
    rownames(dense_matrix2) <- bin_positions
    colnames(dense_matrix2) <- bin_positions
    
    # Populate matrices
    for (i in 1:nrow(contact_matrix1)) {
      bin1 <- (contact_matrix1[i, 1] / resolution) + 1
      bin2 <- (contact_matrix1[i, 2] / resolution) + 1
      count <- contact_matrix1[i, 3]
      if (bin1 <= n_bins && bin2 <= n_bins) {
        dense_matrix1[bin1, bin2] <- count
        dense_matrix1[bin2, bin1] <- count
      }
    }
    
    for (i in 1:nrow(contact_matrix2)) {
      bin1 <- (contact_matrix2[i, 1] / resolution) + 1
      bin2 <- (contact_matrix2[i, 2] / resolution) + 1
      count <- contact_matrix2[i, 3]
      if (bin1 <= n_bins && bin2 <= n_bins) {
        dense_matrix2[bin1, bin2] <- count
        dense_matrix2[bin2, bin1] <- count
      }
    }
    
    # Run differential TAD analysis
    results <- TADCompare(dense_matrix1, dense_matrix2, resolution = resolution, window = window_size)
    
    # Extract and transform data for plotting
    plot_data <- results$Count_Plot$data %>%
      filter(!is.na(Count)) %>%
      mutate(Chromosome = chr, Comparison = paste(tissue1, "vs", tissue2))
    
    # Append to the all_plot_data dataframe
    all_plot_data <- rbind(all_plot_data, plot_data)
  }
}

# Calculate total counts for each chromosome & comparison
chromosome_totals <- all_plot_data %>%
  group_by(Chromosome, Comparison) %>%
  summarise(Total_Count = sum(Count))

# Merge back into the original data
all_plot_data <- all_plot_data %>%
  left_join(chromosome_totals, by = c("Chromosome", "Comparison")) %>%
  mutate(Percentage_Per_Chromosome = (Count / Total_Count) * 100) %>%
  select(-Total_Count)  # Remove the total count column

# Define the order for the 'Type' factor
all_plot_data$Type <- factor(all_plot_data$Type, levels = c("Split", "Complex", "Merge", "Shifted", "Strength Change", "NA", "Non-Differential"))

print(paste("Processing comparison:", tissue1, "vs", tissue2))

# Save the combined data to a file
write.table(all_plot_data, "/Users/users/Desktop/Micro-C/tables_for_paper/diff_tads/tad_changes_eupsc_allchrs_allcomparisons_100k_1mbwindow.txt", row.names = FALSE, quote = FALSE, sep = "\t")

# Ensure each comparison sums to 100% per genome
comparison_totals <- all_plot_data %>%
  group_by(Comparison) %>%
  summarise(Total_Count = sum(Count), .groups = "drop")

# Normalize percentage per comparison (not per chromosome)
all_plot_data <- all_plot_data %>%
  left_join(comparison_totals, by = "Comparison") %>%
  mutate(Percentage_Per_Comparison = (Count / Total_Count) * 100) %>%
  select(-Total_Count)  # Remove the total count column

# Ensure factor ordering for comparisons
all_plot_data$Comparison <- factor(all_plot_data$Comparison, levels = unique(all_plot_data$Comparison))

# Define the correct order for Type
all_plot_data$Type <- factor(all_plot_data$Type, 
                             levels = c("Split", "Complex", "Merge", "Shifted", "Strength Change", "Non-Differential"))

# Define colors for each category
category_colors <- c("Complex" = "#b2df8a", 
                     "Merge" = "#8FCE8F", 
                     "Non-Differential" = "lightgoldenrod", 
                     "Shifted" = "#33a02c", 
                     "Split" = "#e0f8e0", 
                     "Strength Change" =  "#006400", 
                     "NA" = "lightgrey")

# Create the stacked bar plot for all comparisons
p <- ggplot(all_plot_data, aes(x = Comparison, y = Percentage_Per_Comparison, fill = Type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = category_colors) +
  labs(title = "Differential TADs across tissues",
       x = "Tissue Comparison",
       y = "Percentage") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold", margin = margin(b = 15)),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Rotate tissue labels
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    plot.margin = margin(t = 20, r = 12, b = 10, l = 10, unit = "pt")
  ) +
  scale_y_continuous(limits = c(0, 100.9))


# Print the final combined plot
print(p)

ggsave("/Users/users/Desktop/Micro-C/figs_for_paper/Figure4/bar_tad_changes_eupsc_allcomparisons_100k_1mbwindow_unlabeled.tiff", 
       plot = p, width = 10, height = 6, dpi = 200)

#With correct tissue names and ordered----
# Define a mapping for tissue names
tissue_name_map <- c(
  "20" = "stage 20",
  "25" = "stage 25",
  "29" = "stage 29"
)

# Replace tissue names in the Comparison column
all_plot_data$Comparison <- gsub("20", "stage 20", all_plot_data$Comparison)
all_plot_data$Comparison <- gsub("25", "stage 25", all_plot_data$Comparison)
all_plot_data$Comparison <- gsub("29", "stage 29", all_plot_data$Comparison)

# Define custom order for comparisons with full names
comparison_order <- c(
  "stage stage 20 vs stage stage 29", "stage stage 20 vs stage stage 25", "stage stage 29 vs stage stage 25"
)

# Define new names for the comparisons
new_names <- c(
  "stage stage 20 vs stage stage 29" = "Stage 20 vs.\nstage 29",
  "stage stage 20 vs stage stage 25" = "Stage 20 vs.\nstage 25",
  "stage stage 29 vs stage stage 25" = "Stage 25 vs.\nstage 29")

# Ensure 'Comparison' column is a factor and assign the new levels and labels
all_plot_data_fixed <- all_plot_data
all_plot_data_fixed$Comparison <- factor(all_plot_data_fixed$Comparison, levels = comparison_order)
all_plot_data_fixed$Comparison <- recode(all_plot_data_fixed$Comparison, !!!new_names)

# Create the plot
p_fixed <- ggplot(all_plot_data_fixed, aes(x = Comparison, y = Percentage_Per_Comparison, fill = Type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = category_colors) +
  labs(title = "Differential TADs across developmental stages",
       x = "Tissue comparison",
       y = "Percentage") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold", margin = margin(b = 15)),
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Rotate x labels for readability
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    plot.margin = margin(t = 20, r = 12, b = 10, l = 10, unit = "pt")
  ) +
  scale_y_continuous(limits = c(0, 100.9))

# Print the final plot
print(p_fixed)

# Save the final combined plot
ggsave("/Users/users/Desktop/Micro-C/figs_for_paper/Figure4/bar_tad_changes_eupsc_allcomparisons_100k_1mbwindow.tiff", 
       plot = p_fixed, width = 10, height = 6, dpi = 200)

