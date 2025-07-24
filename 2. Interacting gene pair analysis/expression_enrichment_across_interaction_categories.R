#Check expression of different interaction categories

# Clear workspace
rm(list = ls())

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(broom)
library(purrr)
library(grid)
library(VennDiagram)
library(pheatmap)

# Reading the (TPM normalized) gene expression file
expression_data <- read.delim("/Users/users/Desktop/Micro-C/expression_data/lisa_spencer_exp/TPM_countDataMatrixNotFilt.txt", header = TRUE)

# Define the new tissue names in the same order as the sample IDs
tissue_names <- c("gene_id", "B1 arm", "Mantle", "Left gill", "Hectocotylus (A1)", "Skin", 
                  "Testes", "Left optic lobe", "B4 arm", "Right tentacle", "Right gill", "Supraesophogeal lobe", 
                  "Subesophogeal lobe", "Left optic lobe", "Right optic lobe", "Left white body", "Right white body", 
                  "Left gill", "Right gill", "Testes", "Subesophogeal lobe", "Left optic lobe", "Left white body", 
                  "Right white body", "Left gill", "Right gill", "Ovaries", "ANG", "Central core", "Central brain", 
                  "Left optic lobe", "Right optic lobe", "Left white body", "Right white body", "Ovaries", "ANG", 
                  "Central core", "ANG", "Central core", "ANG", "Central core", "ANG")

# Assign the column names to the data frame
colnames(expression_data) <- tissue_names

# Convert all columns except 'gene_id' to numeric and log-transform (adding a small pseudocount to avoid log(0))
expression_data[,-1] <- lapply(expression_data[,-1], function(x) {
  log(as.numeric(as.character(x)) + 0.01)
})

# Replace NA values with 0
expression_data[,-1][is.na(expression_data[,-1])] <- 0

# Convert all columns except 'gene_id' to numeric
expression_data[,-1] <- lapply(expression_data[,-1], as.numeric)

# Convert data to long format
long_data <- expression_data %>%
  pivot_longer(
    cols = -gene_id,
    names_to = "Tissue",
    values_to = "Expression"
  ) %>%
  mutate(Tissue = gsub("\\.\\d+$", "", Tissue)) # Extract base names (remove numeric suffixes)

# Compute average expression per tissue to simplify
averaged_expression_data <- long_data %>%
  group_by(gene_id, Tissue) %>%
  summarise(Average_Expression = mean(Expression, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(
    names_from = Tissue,
    values_from = Average_Expression
  )

# Load interaction data
interaction_data <- read.delim("/Users/users/Desktop/Micro-C/topology_strength_analysis/eupsc_100k/409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof.txt")

# Split the gene pairs into separate columns and isolate gene IDs for the species
interaction_data <- interaction_data %>%
  separate(Orth_pair_based_on_eupsc_interaction_matrix, into = c("Gene1", "Gene2"), sep = ",") %>%
  mutate(Gene1 = sub(";.*", "", Gene1),
         Gene2 = sub(";.*", "", Gene2))

# Merge gene expression data with interaction data
merged_data <- interaction_data %>%
  inner_join(averaged_expression_data, by = c("Gene1" = "gene_id")) %>%
  inner_join(averaged_expression_data, by = c("Gene2" = "gene_id"), suffix = c(".Gene1", ".Gene2"))

#Seperate gene 1 and 2 data and put them on top of eachother----
# Extract columns related to Gene1
gene1_columns <- grep("\\.Gene1$", names(merged_data), value = TRUE)

# Reshape the data for Gene1 and rename Gene1 to Gene
gene1_data <- merged_data %>%
  dplyr::select(Gene1, Interaction_status, all_of(gene1_columns)) %>%
  pivot_longer(
    cols = all_of(gene1_columns),
    names_to = "Tissue",
    values_to = "Expression"
  ) %>%
  mutate(
    Tissue = gsub("\\.Gene1$", "", Tissue),
    Gene = Gene1  # Rename Gene1 to Gene
  ) %>%
  dplyr::select(-Gene1)  # Drop the old Gene1 column

# Extract columns related to Gene2
gene2_columns <- grep("\\.Gene2$", names(merged_data), value = TRUE)

# Reshape the data for Gene2 and rename Gene2 to Gene
gene2_data <- merged_data %>%
  dplyr::select(Gene2, Interaction_status, all_of(gene2_columns)) %>%
  pivot_longer(
    cols = all_of(gene2_columns),
    names_to = "Tissue",
    values_to = "Expression"
  ) %>%
  mutate(
    Tissue = gsub("\\.Gene2$", "", Tissue),
    Gene = Gene2  # Rename Gene2 to Gene
  ) %>%
  dplyr::select(-Gene2)  # Drop the old Gene2 column

# Combine gene1_data and gene2_data
combined_data <- bind_rows(gene1_data, gene2_data)

# Reorder columns to place Gene first
combined_data <- combined_data %>%
  dplyr::select(Gene, everything())  # Move Gene column to the front (for heatmaps later)

# Remove duplicate rows
combined_data <- combined_data %>%
  distinct()

# ---- REMOVE GENES FROM "not_interacting_any_species" IF THEY APPEAR IN OTHER CATEGORIES ----
# Step 1: Identify genes that are in any interacting category
interacting_genes <- combined_data %>%
  filter(Interaction_status != "not_interacting_any_species") %>%
  pull(Gene) %>%
  unique()

# Step 2: Remove genes from "not_interacting_any_species" if they also appear in interacting categories
combined_data <- combined_data %>%
  filter(
    Interaction_status != "not_interacting_any_species" |  # Keep all interacting genes
      (Interaction_status == "not_interacting_any_species" & !Gene %in% interacting_genes)  # Keep only exclusive non-interacting genes
  )


#Normalisation checks----
central_brain_data <- combined_data %>% filter(Tissue == "Central brain") %>% pull(Expression)
testes_data <- combined_data %>% filter(Tissue == "Testes") %>% pull(Expression)
ang_data <- combined_data %>% filter(Tissue == "ANG") %>% pull(Expression)

# Plotting the histogram in R
hist(central_brain_data, 
     main = "Histogram of gene expression in central brain tissue", 
     xlab = "Expression level", 
     col = "blue", 
     border = "black", 
     breaks = 30)


# Plotting the histogram in R
hist(testes_data, 
     main = "Histogram of gene expression in testes tissue", 
     xlab = "Expression level", 
     col = "blue", 
     border = "black", 
     breaks = 30)


# Plotting the histogram in R
hist(ang_data, 
     main = "Histogram of gene expression in ANG tissue", 
     xlab = "Expression level", 
     col = "blue", 
     border = "black", 
     breaks = 30)

# Enrichment analysis----
# Perform pairwise Wilcoxon tests within each tissue comparing Interaction_status categories
pairwise_wilcox_results <- combined_data %>%
  group_by(Tissue) %>%
  group_split() %>%
  map_dfr(~ {
    data_subset <- .
    
    # Get unique Interaction_status categories
    status_levels <- unique(data_subset$Interaction_status)
    
    # Generate unique pairs (only one direction)
    status_pairs <- combn(status_levels, 2, simplify = FALSE)
    
    # Perform pairwise Wilcoxon tests
    results <- map_dfr(status_pairs, function(pair) {
      status1 <- pair[1]
      status2 <- pair[2]
      
      subset1 <- data_subset %>% filter(Interaction_status == status1) %>% pull(Expression)
      subset2 <- data_subset %>% filter(Interaction_status == status2) %>% pull(Expression)
      
      wilcox_test <- wilcox.test(subset1, subset2, exact = FALSE)
      
      tibble(
        Interaction_status1 = status1,
        Interaction_status2 = status2,
        p_value = wilcox_test$p.value,
        Tissue = unique(data_subset$Tissue)
      )
    })
    
    return(results)
  })

# Benjamini-Hochberg FDR
pairwise_wilcox_results <- pairwise_wilcox_results %>%
  mutate(p_adjusted = p.adjust(p_value, method = "BH"))  # Benjamini-Hochberg FDR


# Save results
write.table(pairwise_wilcox_results, file = "/Users/users/Desktop/Micro-C/tables_for_paper/expression_of_spatiosyntenies/tissue_expression_wilcox_eupsc_prioritise_both_unique_neither_with_sof.txt", row.names = FALSE, quote = FALSE, sep = "\t")


# Plot boxplots for each tissue, with each Interaction_status represented as a boxplot within each tissue----

# Define the custom colors for filling
interaction_cols <- c(
  "interacting_all_species" = adjustcolor("#FF7878", alpha.f = 0.8), 
  "interacting_deca_only" = adjustcolor("palegreen", alpha.f = 0.8), 
  "interacting_octbi_only" = adjustcolor("mediumpurple1", alpha.f = 0.8), 
  "not_interacting_any_species" = adjustcolor("orange", alpha.f = 0.8)
)

# Order
combined_data$Interaction_status <- factor(combined_data$Interaction_status, 
                                           levels = c("interacting_all_species", 
                                                      "interacting_deca_only", 
                                                      "interacting_octbi_only", 
                                                      "not_interacting_any_species"))

# Define the desired order of tissues
desired_order <- c("Subesophogeal lobe", "Supraesophogeal lobe", "Central brain",
                   "Left optic lobe", "Right optic lobe", "Hectocotylus (A1)",
                   "B1 arm", "B4 arm", "Right tentacle", "Left white body",
                   "Right white body", "ANG", "Central core", "Skin", "Mantle",
                   "Left gill", "Right gill", "Ovaries", "Testes")
# Create the plot
tissue_boxplots <- combined_data %>%
  ggplot(aes(x = Interaction_status, y = Expression, fill = Interaction_status)) +
  geom_boxplot(color = "black",  # Outline color for the boxplots
               outlier.shape = 19,  # Shape of outliers
               outlier.size = 0.5,  # Size of outliers
               outlier.colour = "black") +  # Color of outliers
  facet_wrap(~ factor(Tissue, levels = desired_order), scales = "free_y") +
  scale_fill_manual(values = interaction_cols) +  # Apply custom fill colors
  theme_bw() +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis labels
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.ticks.x = element_blank(),  # Remove x-axis tick marks
    panel.border = element_rect(color = "black", fill = NA),  # Black border for panels
    legend.position = "right",  # Position legend on the right
    strip.text = element_text(size = 6),  # Smaller text for facet labels
    plot.title = element_text(size = 10)  # Smaller text for plot title
  ) +
  labs(
    title = "E. scolopes gene expression per tissue for each category of interaction status",
    x = NULL,  # Set x-axis label to NULL
    y = "Expression (logTPM)"
  )

# Print the plot
print(tissue_boxplots)

# Save the plot
ggsave(tissue_boxplots, file = "/Users/users/Desktop/Micro-C/figs_for_paper/Supplementary/tissue_expression/tissue_expression_boxplot_prioritise_both_unique_neither_eupsc_with_sof.tiff", height = 6, width = 9)


#Heatmaps----

# Get unique interaction statuses
interaction_statuses <- unique(combined_data$Interaction_status)

#NOT CLUSTERED BY GENE----

# Define a function to create and plot heatmaps.---
plot_heatmap_for_status <- function(data, status) {
  # Subset data for the current interaction status and drop the Interaction_status column
  subset_data <- data %>%
    filter(Interaction_status == status) %>%
    dplyr::select(-Interaction_status)
  
  # Reshape data into matrix format
  heatmap_matrix <- subset_data %>%
    pivot_wider(
      names_from = Tissue,
      values_from = Expression
    )
  
  # Convert the reshaped dataframe to a matrix
  heatmap_matrix <- as.data.frame(heatmap_matrix)
  row.names(heatmap_matrix) <- heatmap_matrix$Gene
  heatmap_matrix <- heatmap_matrix %>%
    dplyr::select(-Gene) %>%
    as.matrix()
  
  # Ensure no NA values for heatmap plotting
  heatmap_matrix[is.na(heatmap_matrix)] <- 0
  
  # Create heatmap with default colors, dendrograms for columns, and no grey outlines
  pheatmap(
    heatmap_matrix,
    cluster_rows = FALSE,   # Do not cluster rows (genes)
    cluster_cols = TRUE,    # Cluster columns (tissues)
    show_rownames = FALSE,  # Do not show row names (genes)
    show_colnames = TRUE,
    main = paste("Heatmap for", status), # Add title for each heatmap
    treeheight_row = 0,    # Hide row dendrogram while keeping clustering
    treeheight_col = 50,    # Adjust height of the column dendrogram
    border_color = NA,      # Remove borders around cells
    legend = TRUE           # Show legend
  )
}

# Plot heatmaps for each interaction status
for (status in interaction_statuses) {
  plot_heatmap_for_status(combined_data, status)
}


#CLUSTERED BY GENE - you must remove rows that have zero variance to do this.----

plot_heatmap_for_status_gene_cluster <- function(data, status, output_dir = "/Users/users/Desktop/Micro-C/figs_for_paper/Supplementary/tissue_expression/") {
  # Create output directory if it does not exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Subset data for the current interaction status and drop the Interaction_status column
  subset_data <- data %>%
    filter(Interaction_status == status) %>%
    dplyr::select(-Interaction_status)
  
  # Reshape data into matrix format
  heatmap_matrix <- subset_data %>%
    pivot_wider(
      names_from = Tissue,
      values_from = Expression
    )
  
  # Convert the reshaped dataframe to a matrix
  heatmap_matrix <- as.data.frame(heatmap_matrix)
  row.names(heatmap_matrix) <- heatmap_matrix$Gene
  heatmap_matrix <- heatmap_matrix %>%
    dplyr::select(-Gene) %>%
    as.matrix()
  
  # Ensure no NA values for heatmap plotting
  heatmap_matrix[is.na(heatmap_matrix)] <- 0
  
  # Debugging: Check if there are any non-numeric values
  if (any(!is.numeric(heatmap_matrix))) {
    print("Non-numeric data found in matrix.")
  }
  
  # Debugging: Check for zero variance rows
  zero_var_rows <- apply(heatmap_matrix, 1, var) == 0
  if (any(zero_var_rows)) {
    print("Zero variance rows found.")
    print(row.names(heatmap_matrix)[zero_var_rows])
  }
  
  # Remove zero variance rows if necessary
  heatmap_matrix <- heatmap_matrix[!zero_var_rows, ]
  
  # Define the filename for the heatmap
  filename <- file.path(output_dir, paste("tissue_expression_heatmap_eupsc_prioritise_both_unique_neither_", status, "_with_sof.tiff", sep = ""))
  
  # Create heatmap and save to file
  heatmap_plot <- pheatmap(
    heatmap_matrix,
    cluster_rows = TRUE,   # Perform clustering on rows
    cluster_cols = TRUE,    # Cluster columns (tissues)
    show_rownames = FALSE,  # Do not show row names (genes)
    show_colnames = TRUE,     
    main = paste("Heatmap for", status), # Add title for each heatmap
    treeheight_row = 0,    # Adjust height of the row dendrogram (0 to remove it)
    treeheight_col = 50,    # Adjust height of the column dendrogram
    border_color = NA,      # Remove borders around cells
    legend = TRUE,          # Show legend
    filename = filename,     # Save to file
    breaks = seq(-8, 8, length.out = 101),# Force axes
  )
  
  # Print the heatmap to the console
  grid::grid.draw(heatmap_plot$gtable)
}

# Plot and save heatmaps for each interaction status
for (status in interaction_statuses) {
  plot_heatmap_for_status_gene_cluster(combined_data, status)
}


#Venn diagrams of each category ----

# Close the current plot device if open
if (dev.cur() > 1) dev.off()

# Function to get gene lists by interaction status
get_gene_lists <- function(data) {
  interaction_statuses <- unique(data$Interaction_status)
  gene_lists <- list()
  
  for (status in interaction_statuses) {
    gene_lists[[status]] <- data %>%
      filter(Interaction_status == status) %>%
      pull(Gene) %>%
      unique()
  }
  
  return(gene_lists)
}

# Get the gene lists
gene_lists <- get_gene_lists(combined_data)

# Reorder gene_lists so it follows the order of interaction_cols
gene_lists <- gene_lists[names(interaction_cols)]

# Plot the Venn diagram
tissue_venn <- venn.diagram(
  x = gene_lists,
  category.names = names(gene_lists),
  filename = NULL,    # Display the plot in R console instead of saving to a file
  output = TRUE,
  main = NULL,
  col = "black",
  fill = c(interaction_cols), # Adjust colors as needed
  alpha = 0.8,        # Transparency of colors
  cex = 1.2,            # Larger font size for counts
  cat.cex = 0,        # Hide category labels
  cat.col = c(interaction_cols), # Adjust colors as needed
  fontfamily = "sans")

# Display the Venn diagram
grid.draw(tissue_venn)

# Save the Venn diagram to a file
ggsave(tissue_venn, file = "/Users/users/Desktop/Micro-C/figs_for_paper/Supplementary/tissue_expression/tissue_expression_venn_eupsc_prioritise_both_unique_neither_with_sof.tiff", width = 8, height = 6.5)

#Tau calculation
# Compute Tau for each gene per interaction category
calculate_tau <- function(expression_values) {
  x_max <- max(expression_values, na.rm = TRUE)
  if (x_max <= 0) return(0)  # Avoid division by zero or negative max
  n <- sum(!is.na(expression_values))  # Count only non-NA values
  if (n <= 1) return(0)  # If only one tissue, Tau is 0
  
  tau <- sum(1 - (expression_values / x_max), na.rm = TRUE) / (n - 1)
  return(max(0, min(tau, 1)))  # Ensure Tau is between 0 and 1
}

# Compute Tau for each gene in each interaction category
tau_results <- combined_data %>%
  group_by(Gene, Interaction_status) %>%
  summarise(Tau = calculate_tau(Expression), .groups = 'drop')

# Compute summary statistics per interaction category
summary_tau <- tau_results %>%
  group_by(Interaction_status) %>%
  summarise(
    Mean_Tau = mean(Tau, na.rm = TRUE),
    Median_Tau = median(Tau, na.rm = TRUE),
    .groups = 'drop'
  )

# Perform Wilcoxon pairwise comparisons between interaction categories
pairwise_wilcox <- pairwise.wilcox.test(tau_results$Tau, tau_results$Interaction_status, p.adjust.method = "BH")
pairwise_wilcox_unadjusted <- pairwise.wilcox.test(tau_results$Tau, tau_results$Interaction_status, p.adjust.method = "none")

# Convert matrix to long format while ensuring all comparisons are kept
wilcox_results <- as.data.frame(as.table(pairwise_wilcox$p.value))
colnames(wilcox_results) <- c("Comparison", "Interaction_status", "p_adj")

wilcox_results_unadjusted <- as.data.frame(as.table(pairwise_wilcox_unadjusted$p.value))
colnames(wilcox_results_unadjusted) <- c("Comparison", "Interaction_status", "p")

# Merge p-values and adjusted p-values
wilcox_results <- wilcox_results_unadjusted %>% 
  left_join(wilcox_results, by = c("Comparison", "Interaction_status")) %>%
  drop_na()  # Remove NA values

# Save results to files
write.table(tau_results, 
            file = "/Users/users/Desktop/Micro-C/tables_for_paper/expression_of_spatiosyntenies/tau_tissue_specificity.txt", 
            row.names = FALSE, 
            quote = FALSE, 
            sep = "\t")

write.table(summary_tau, 
            file = "/Users/users/Desktop/Micro-C/tables_for_paper/expression_of_spatiosyntenies/summary_tau_per_category.txt", 
            row.names = FALSE, 
            quote = FALSE, 
            sep = "\t")

write.table(wilcox_results, 
            file = "/Users/users/Desktop/Micro-C/tables_for_paper/expression_of_spatiosyntenies/wilcox_tau_comparisons.txt", 
            row.names = FALSE, 
            quote = FALSE, 
            sep = "\t")

# Print results
print(head(tau_results))
print(summary_tau)
print(wilcox_results)


#Tau for development
# Clear environment
rm(list = ls())

# Load libraries
library(dplyr)
library(tidyr)

# Load developmental expression data
dev_expr_data <- read.delim("/Users/users/Desktop/Micro-C/expression_data/Hannah_TPM_normalizedcounts_development.txt", header = TRUE)

# Log-transform with pseudocount
dev_expr_data[,-1] <- lapply(dev_expr_data[,-1], function(x) log(as.numeric(as.character(x)) + 0.01))

# Replace NA with 0 (optional)
dev_expr_data[,-1][is.na(dev_expr_data[,-1])] <- 0

# Load interaction data with Interaction_status and Gene columns
interaction_data <- read.delim("/Users/users/Desktop/Micro-C/tables_for_paper/expression_of_spatiosyntenies/tau_tissue_specificity.txt")

# Make sure only needed columns are kept
gene_status <- interaction_data %>%
  dplyr::select(Gene, Interaction_status) %>%
  distinct()

# Filter expression data to matching genes
dev_expr_data_filtered <- dev_expr_data %>%
  filter(gene_id %in% gene_status$Gene)

# Pivot to long format and join with interaction status
dev_long <- dev_expr_data_filtered %>%
  pivot_longer(-gene_id, names_to = "Stage", values_to = "Expression") %>%
  inner_join(gene_status, by = c("gene_id" = "Gene"))

# Tau function
calculate_tau <- function(expression_values) {
  x_max <- max(expression_values, na.rm = TRUE)
  if (x_max <= 0) return(0)
  n <- sum(!is.na(expression_values))
  if (n <= 1) return(0)
  tau <- sum(1 - (expression_values / x_max), na.rm = TRUE) / (n - 1)
  return(max(0, min(tau, 1)))
}

# Compute Tau per gene per interaction category
tau_results_dev <- dev_long %>%
  group_by(gene_id, Interaction_status) %>%
  summarise(Tau = calculate_tau(Expression), .groups = 'drop')

# Summary statistics per interaction status
summary_tau_dev <- tau_results_dev %>%
  group_by(Interaction_status) %>%
  summarise(
    Mean_Tau = mean(Tau, na.rm = TRUE),
    Median_Tau = median(Tau, na.rm = TRUE),
    .groups = 'drop'
  )

# Save results
write.table(tau_results_dev,
            file = "/Users/users/Desktop/Micro-C/tables_for_paper/expression_of_spatiosyntenies/tau_development_by_interaction_status.txt",
            row.names = FALSE, quote = FALSE, sep = "\t")

write.table(summary_tau_dev,
            file = "/Users/users/Desktop/Micro-C/tables_for_paper/expression_of_spatiosyntenies/summary_tau_development_by_interaction_status.txt",
            row.names = FALSE, quote = FALSE, sep = "\t")

# Print
print(summary_tau_dev)