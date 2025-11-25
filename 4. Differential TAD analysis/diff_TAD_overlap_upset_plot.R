#Note when a TAD appears in only some files, it's considered differential, even if the TADs are non-differential in the files it does appear in.
#E.g.  if a TAD appears in three out of six files and is marked as "Non-Differential" in those three files, but does not appear in the other three files, it will still be considered differential and included in the plot.
#Any TADs with NAs in the differential column will be removed from the files
rm(list = ls())

# Load required libraries
library(UpSetR)
library(dplyr)
library(tidyr)
library(ggplot2)

# List all comparison files
file_paths <- Sys.glob("/Users/users/Desktop/Micro-C/tables_for_paper/diff_tads/differential_tads_eupsc*")

# Read all files and extract differential TAD boundaries
tad_list <- lapply(file_paths, function(file) {
  df <- read.delim(file, header = TRUE, sep = "\t", colClasses = c("Boundary" = "character"))
  df$Comparison <- gsub("differential_tads_|_allchrs_100k_1mbwindow.txt", "", basename(file))  # Extract comparison name
  
  # Create TAD_ID before filtering
  df <- df %>%
    mutate(TAD_ID = paste(Chromosome, Boundary, sep = "_"))  # Create unique TAD ID
  
  return(df)
})

# Combine all TAD data into one table
all_tads <- bind_rows(tad_list)

# Remove TADs with "Non-Differential" in **all files** and remove TADs with NA in any file
filtered_tads <- all_tads %>%
  group_by(TAD_ID) %>%
  filter(!all(Differential == "Non-Differential")) %>%  # Keep if not "Non-Differential" in all files
  filter(!any(is.na(Differential))) %>%  # Remove TADs with NA in any comparison
  ungroup()

# Ensure each TAD_ID appears once per comparison
filtered_tads <- filtered_tads %>%
  mutate(TAD_ID = paste(Chromosome, Boundary, sep = "_")) %>%  # Create unique TAD ID again
  distinct(TAD_ID, Comparison)

# Convert to wide format: Each row is a TAD, columns are comparisons (1 = present, 0 = absent)
tad_matrix <- filtered_tads %>%
  pivot_wider(names_from = Comparison, values_from = Comparison, values_fill = NA) %>%
  mutate(across(-TAD_ID, ~ ifelse(is.na(.), 0, 1))) %>%  # Convert presence/absence to 0/1
  as.data.frame()  # Ensure it's a standard data frame, not a tibble

# Convert all comparison columns to integer (UpSetR requires numeric input)
tad_matrix[, -1] <- lapply(tad_matrix[, -1], as.integer)

# Generate the UpSet plot
upset(
  tad_matrix[, -1],  # Remove TAD_ID and use only binary matrix
  sets = colnames(tad_matrix)[-1], 
  mainbar.y.label = "Differential TAD overlap",
  sets.x.label = "Total differential TADs per comparison",
  order.by = "freq",
  keep.order = TRUE
)

# Ordered and coloured ----

# Calculate total count per comparison
set_sizes <- colSums(tad_matrix[, -1])  # Sum each column except TAD_ID

# Order sets from smallest to largest for the left barplot
ordered_sets <- names(sort(set_sizes, decreasing = FALSE))  # Now in ascending order

# Reorder `tad_matrix` columns based on ascending order of set sizes
tad_matrix <- tad_matrix[, c("TAD_ID", ordered_sets)]

# Convert to numeric to ensure UpSetR compatibility
tad_matrix[, -1] <- lapply(tad_matrix[, -1], as.integer)

# Compute total TAD counts per comparison (left bars)
set_sizes <- colSums(tad_matrix[, -1])  # Exclude TAD_ID column

# Get the number of left bars (tissue comparisons)
num_left_bars <- length(set_sizes)

# Generate a purple color gradient for the left bars
left_bar_colors <- colorRampPalette(c("#de425b", "#82204a"))(num_left_bars)
#"#CCDDCC", "#88CCAA", "#44AA99", "#228866", "#114433"
# "#88CCAA", "#117733", "#004422"

# Apply the UpSet plot
upset_plot <- upset(
  tad_matrix[, -1],  # Exclude TAD_ID column
  sets = names(set_sizes),  # Use the tissue comparisons (left bars)
  mainbar.y.label = "Differential TAD overlap",
  sets.x.label = "Total differential TADs per comparison",
  order.by = "freq",  # Keep top bars in descending order
  keep.order = TRUE,
  sets.bar.color = left_bar_colors,  # Apply purple gradient to left bars
  main.bar.color = "#33a02c",  # Set the top bars to grey
  text.scale = 1.5  # Increase text size for labels and titles
)

# Print the plot
print(upset_plot)

#Renamed axes
new_colnames <- c(
  "eupsc_20_29_allchrs_1mbwindow.txt" = "Stage 20 vs. stage 29",
  "eupsc_25_29_allchrs_1mbwindow.txt" = "Stage 25 vs. stage 29",
  "eupsc_20_25_allchrs_1mbwindow.txt" = "Stage 20 vs. stage 25"
)

# Rename the columns in tad_matrix
colnames(tad_matrix) <- ifelse(
  colnames(tad_matrix) %in% names(new_colnames),
  new_colnames[colnames(tad_matrix)],
  colnames(tad_matrix)
)

upset_plot_final <- upset(
  tad_matrix[, -1],  # Exclude TAD_ID column
  sets = colnames(tad_matrix[, -1]),  # Use the renamed column names directly
  mainbar.y.label = "Differential TAD overlap",
  sets.x.label = "Total differential TADs per comparison",
  order.by = "freq",  # Keep top bars in descending order
  keep.order = TRUE,
  sets.bar.color = left_bar_colors,  # Apply purple gradient to left bars
  main.bar.color = "#005500",  # Set the top bars to grey
  text.scale = 1.49  # Ensure the text scale is consistent
)


# Save as TIFF
tiff("/Users/users/Desktop/Micro-C/figs_for_paper/Figure4/diff_tad_tissue_overlap_eupsc_100k_1mbwindow.tiff", width = 14.1, height = 6, units = "in", res = 500)
print(upset_plot_final)  # This will plot to the TIFF device
dev.off()


#Count on your terminal he number of differential TADS

#grep -v NA /Users/users/Desktop/Micro-C/tables_for_paper/diff_tads/differential_tads_eupsc_cb_ar_allchrs_1mbwindow.txt | grep 'Non-Differential' | wc -l
#   1921

#grep -v NA /Users/users/Desktop/Micro-C/tables_for_paper/diff_tads/differential_tads_eupsc_cb_ar_allchrs_1mbwindow.txt | grep -v 'Non-Differential' | wc -l
#     730

#grep -v NA /Users/users/Desktop/Micro-C/tables_for_paper/diff_tads/differential_tads_eupsc_*_*_allchrs_1mbwindow.txt | grep 'Non-Differential' | wc -l
#  11411

#grep -v NA /Users/users/Desktop/Micro-C/tables_for_paper/diff_tads/differential_tads_eupsc_*_*_allchrs_1mbwindow.txt | grep -v 'Non-Differential' | wc -l
#    4703
