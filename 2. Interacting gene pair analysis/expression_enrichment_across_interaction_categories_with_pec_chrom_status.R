# # Check expression of different interaction and P. maximus chromosome status categories by plotting boxplots, heatmaps, Venn diagrams, and performing significance tests and calculation of tau (tissue specificity) For E. scolopes tissues

rm(list = ls())

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(broom)
library(purrr)
library(grid)
library(pheatmap)
library(tibble)

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

# Replace NA values with 0 (if there are any - there aren't but you should really do this before...)
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

# Split the gene pairs into separate columns and isolate IDs for the species
interaction_data <- interaction_data %>%
  separate(Orth_pair_based_on_eupsc_interaction_matrix, into = c("Gene1", "Gene2"), sep = ",") %>%
  mutate(Gene1 = sub(";.*", "", Gene1),
         Gene2 = sub(";.*", "", Gene2))

# Merge gene expression data with interaction data
merged_data <- interaction_data %>%
  inner_join(averaged_expression_data, by = c("Gene1" = "gene_id")) %>%
  inner_join(averaged_expression_data, by = c("Gene2" = "gene_id"), suffix = c(".Gene1", ".Gene2"))

# Create a Combined_status column based on Pecten_chromosome_status and Interaction_status
merged_data <- merged_data %>%
  mutate(Combined_status = paste(Interaction_status, Pecten_chromosome_status, sep = "_"))

# Extract columns related to Gene1 and Gene2
gene1_columns <- grep("\\.Gene1$", names(merged_data), value = TRUE)
gene2_columns <- grep("\\.Gene2$", names(merged_data), value = TRUE)

# Reshape the data for Gene1 and Gene2 and combine them
combined_data <- bind_rows(
  merged_data %>%
    dplyr::select(Gene1, Interaction_status, Combined_status, all_of(gene1_columns)) %>%
    pivot_longer(
      cols = all_of(gene1_columns),
      names_to = "Tissue",
      values_to = "Expression"
    ) %>%
    mutate(Tissue = gsub("\\.Gene1$", "", Tissue), Gene = Gene1) %>%
    dplyr::select(-Gene1),
  
  merged_data %>%
    dplyr::select(Gene2, Interaction_status, Combined_status, all_of(gene2_columns)) %>%
    pivot_longer(
      cols = all_of(gene2_columns),
      names_to = "Tissue",
      values_to = "Expression"
    ) %>%
    mutate(Tissue = gsub("\\.Gene2$", "", Tissue), Gene = Gene2) %>%
    dplyr::select(-Gene2)
)

# Remove duplicate rows
combined_data <- combined_data %>%
  distinct()


#Keep all genes in interacting both, but only keep genes that in interacting neither species if they are ONLY in this category i.e. unique to this category, NOT IN THE OTHER 3! 
#I.e. prioritise interacting both species, and keep interacting neither super unique.

# Filter combined_data to include only:
# - Genes in the "interacting_all_species" category
# - Genes in the "not_interacting_any_species" category that are not in "interacting_all_species" OR interacting either species only


# Step 1: Identify overlapping genes in categories other than "not_interacting_any_species"
interacting_genes <- combined_data %>%
  filter(Interaction_status != "not_interacting_any_species") %>%
  pull(Gene)

# Step 2: Filter the data based on the specified conditions
combined_data <- combined_data %>%
  filter(
    Interaction_status == "interacting_all_species" | # Keep all genes in interacting_all_species
      (Interaction_status == "not_interacting_any_species" & !Gene %in% interacting_genes) # Keep genes in not_interacting_any_species only if they don't overlap with other categories
  )

# Print unique values of Combined_status to check the filtering worked
unique_combined_statuses <- unique(combined_data$Combined_status)
print(unique_combined_statuses)

# Count the number of values in each category of Combined_status to check the filtering worked
combined_status_counts <- combined_data %>%
  group_by(Combined_status) %>%
  summarise(Count = n())

# Print the counts
print(combined_status_counts)


# Enrichment analysis
# Perform pairwise Wilcoxon tests within each tissue for 'Combined_status' categories
pairwise_wilcox_results <- combined_data %>%
  group_by(Tissue) %>%
  group_split() %>%
  map_dfr(~ {
    data_subset <- .
    
    # Get unique Combined_status categories
    status_levels <- unique(data_subset$Combined_status)
    
    # Generate unique pairs (only one direction)
    status_pairs <- combn(status_levels, 2, simplify = FALSE)
    
    # Perform pairwise Wilcoxon tests within each tissue
    results <- map_dfr(status_pairs, function(pair) {
      status1 <- pair[1]
      status2 <- pair[2]
      
      # Subset data for the two Combined_status levels
      subset1 <- data_subset %>% filter(Combined_status == status1) %>% pull(Expression)
      subset2 <- data_subset %>% filter(Combined_status == status2) %>% pull(Expression)
      
      # Perform the Wilcoxon test
      wilcox_test <- wilcox.test(subset1, subset2, exact = FALSE)
      
      # Create a tidy result
      tibble(
        Combined_status1 = status1,
        Combined_status2 = status2,
        p_value = wilcox_test$p.value,
        Tissue = unique(data_subset$Tissue)
      )
    })
    
    return(results)
  })


# Benjamini-Hochberg FDR
pairwise_wilcox_results <- pairwise_wilcox_results %>%
  mutate(p_adjusted = p.adjust(p_value, method = "BH"))

# Save the results to a file
write.table(pairwise_wilcox_results, 
            file = "/Users/users/Desktop/Micro-C/tables_for_paper/expression_of_spatiosyntenies/tissue_expression_wilcox_prioritise_both_unique_neither_combined_status_eupsc_with_sof.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")


# Define colors for each Combined_status
interaction_cols <- c(
  "not_interacting_any_species_diff_pec_chrs" = "lightblue",
  "not_interacting_any_species_same_pec_chrs" = "indianred",
  "interacting_all_species_diff_pec_chrs" = "#377EB8",
  "interacting_all_species_same_pec_chrs" = "#E41A1C"
)

# Define the desired order of each box
status_order <- c("interacting_all_species_same_pec_chrs", 
                  "interacting_all_species_diff_pec_chrs", 
                  "not_interacting_any_species_same_pec_chrs", 
                  "not_interacting_any_species_diff_pec_chrs")


# Convert the Combined_status column to a factor with the desired order
combined_data <- combined_data %>%
  mutate(Combined_status = factor(Combined_status, levels = status_order))

# Define the desired order of tissues
desired_order <- c("Subesophogeal lobe", "Supraesophogeal lobe", "Central brain",
                   "Left optic lobe", "Right optic lobe", "Hectocotylus (A1)",
                   "B1 arm", "B4 arm", "Right tentacle", "Left white body",
                   "Right white body", "ANG", "Central core", "Skin", "Mantle",
                   "Left gill", "Right gill", "Ovaries", "Testes")

# Plot boxplots for each tissue, grouped by Combined_status
tissue_boxplots <- combined_data %>%
  ggplot(aes(x = Combined_status, y = Expression, fill = Combined_status)) +
  geom_boxplot(color = "black",  # Outline color for the boxplots
               outlier.shape = 19,  # Shape of outliers
               outlier.size = 0.5,  # Size of outliers
               outlier.colour = "black") +  # Color of outliers
  facet_wrap(~ factor(Tissue, levels = desired_order), scales = "free_y") +
  theme_bw() +
  scale_fill_manual(values = interaction_cols) +
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
    title = "E. scolopes gene expression per tissue for each category of combined status",
    x = NULL,  # Set x-axis label to NULL
    y = "Expression (logTPM)"
  )

# Print the plot
print(tissue_boxplots)

# Save the plot
ggsave(tissue_boxplots, file = "/Users/users/Desktop/Micro-C/figs_for_paper/Supplementary/tissue_expression_boxplot_prioritise_both_unique_neither_combined_status_eupsc._with_sof.tiff", height = 6, width = 10)



#Heatmaps----

#NOT CLUSTERED BY GENE----

# Define a function to create and plot heatmaps.
plot_heatmap_for_status <- function(data, status) {
  # Subset data for the current interaction status and drop the Combined_status and Interaction_status column or you'll have problems downstream.
  subset_data <- data %>%
    filter(Combined_status == status) %>%
    dplyr::select(-Combined_status, -Interaction_status)
  
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
  
  # Create heatmap without clustering rows
  pheatmap(
    heatmap_matrix,
    cluster_rows = FALSE,   # Do not cluster rows (genes)
    cluster_cols = TRUE,    # Cluster columns (tissues)
    show_rownames = FALSE,  # Do not show row names (genes)
    show_colnames = TRUE, #Don't cluster by row~~~~
    border_color = NA,      # Remove borders around cells
    main = paste("Heatmap for", status) # Add title for each heatmap
  )
}

# Plot heatmaps for each interaction status
for (status in unique_combined_statuses) {
  plot_heatmap_for_status(combined_data, status)
}




#CLUSTERED BY GENE - you must remove rows that have zero variance to do this. ----

plot_heatmap_for_status_gene_cluster <- function(data, status, output_dir = "/Users/users/Desktop/Micro-C/figs_for_paper/Supplementary/tissue_expression") {
  # Create output directory if it does not exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Subset data for the current interaction status and drop both the Combined_status and Interaction_status column or you'll have problems downstream
  subset_data <- data %>%
    filter(Combined_status == status) %>%
    dplyr::select(-Combined_status, -Interaction_status)
  
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
  filename <- file.path(output_dir, paste("tissue_expression_heatmap_prioritise_interacting_both_unique_non_interacting_eupsc_", status, "._with_sof.tiff", sep = ""))
  
  # Create heatmap and save to file
  heatmap_plot <- pheatmap(
    heatmap_matrix,
    cluster_rows = TRUE,   # Perform clustering on rows
    cluster_cols = TRUE,    # Cluster columns (tissues)
    show_rownames = FALSE,  # Do not show row names (genes)
    show_colnames = TRUE,
    main = paste("Heatmap for", status), # Add title for each heatmap
    treeheight_row = 0,    # Hide row dendrogram while keeping clustering
    treeheight_col = 50,    # Adjust height of the column dendrogram
    border_color = NA,      # Remove borders around cells
    legend = TRUE,          # Show legend
    filename = filename,
    breaks = seq(-8, 8, length.out = 101),# Force axes
  )
  
  # Print the heatmap to the console
  grid::grid.draw(heatmap_plot$gtable)
}

# Plot and save heatmaps for each interaction status
for (status in unique_combined_statuses) {
  plot_heatmap_for_status_gene_cluster(combined_data, status)
}

#Venn diagrams of each category ----

# Close the current plot device if open
if (dev.cur() > 1) dev.off()

# Function to get gene lists by interaction status
get_gene_lists <- function(data) {
  combined_statuses <- unique(data$Combined_status)
  gene_lists <- list()
  
  for (status in combined_statuses) {
    gene_lists[[status]] <- data %>%
      filter(Combined_status == status) %>%
      pull(Gene) %>%
      unique()
  }
  
  return(gene_lists)
}

# Get the gene lists
gene_lists <- get_gene_lists(combined_data)

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
ggsave(tissue_venn, file = "/Users/users/Desktop/Micro-C/figs_for_paper/Supplementary/tissue_expression/tissue_expression_venn_eupsc_combined_status._with_sof.tiff", width = 8, height = 6.5)


# Compute Tau for each gene per Combined_status category
calculate_tau <- function(expression_values) {
  x_max <- max(expression_values, na.rm = TRUE)
  if (x_max <= 0) return(0)  # Avoid division by zero or negative max
  n <- sum(!is.na(expression_values))  # Count only non-NA values
  if (n <= 1) return(0)  # If only one tissue, Tau is 0
  
  tau <- sum(1 - (expression_values / x_max), na.rm = TRUE) / (n - 1)
  return(max(0, min(tau, 1)))  # Ensure Tau is between 0 and 1
}

# Compute Tau for each gene in each Combined_status category
tau_results_combined <- combined_data %>%
  group_by(Gene, Combined_status) %>%
  summarise(Tau = calculate_tau(Expression), .groups = 'drop')

# Compute summary statistics per Combined_status category
summary_tau_combined <- tau_results_combined %>%
  group_by(Combined_status) %>%
  summarise(
    Mean_Tau = mean(Tau, na.rm = TRUE),
    Median_Tau = median(Tau, na.rm = TRUE),
    .groups = 'drop'
  )

# Perform Wilcoxon pairwise comparisons between Combined_status categories
pairwise_wilcox_combined <- pairwise.wilcox.test(tau_results_combined$Tau, tau_results_combined$Combined_status, p.adjust.method = "BH")
pairwise_wilcox_unadjusted_combined <- pairwise.wilcox.test(tau_results_combined$Tau, tau_results_combined$Combined_status, p.adjust.method = "none")

# Convert matrix to long format while ensuring all comparisons are kept
wilcox_results_combined <- as.data.frame(as.table(pairwise_wilcox_combined$p.value))
colnames(wilcox_results_combined) <- c("Comparison", "Combined_status", "p_adj")

wilcox_results_unadjusted_combined <- as.data.frame(as.table(pairwise_wilcox_unadjusted_combined$p.value))
colnames(wilcox_results_unadjusted_combined) <- c("Comparison", "Combined_status", "p")

# Merge p-values and adjusted p-values
wilcox_results_combined <- wilcox_results_unadjusted_combined %>% 
  left_join(wilcox_results_combined, by = c("Comparison", "Combined_status")) %>%
  drop_na()  # Remove NA values

# Save results to files
write.table(tau_results_combined, 
            file = "/Users/users/Desktop/Micro-C/tables_for_paper/expression_of_spatiosyntenies/tau_tissue_specificity_combined_status.txt", 
            row.names = FALSE, 
            quote = FALSE, 
            sep = "\t")

write.table(summary_tau_combined, 
            file = "/Users/users/Desktop/Micro-C/tables_for_paper/expression_of_spatiosyntenies/summary_tau_per_combined_status.txt", 
            row.names = FALSE, 
            quote = FALSE, 
            sep = "\t")

write.table(wilcox_results_combined, 
            file = "/Users/users/Desktop/Micro-C/tables_for_paper/expression_of_spatiosyntenies/wilcox_tau_comparisons_combined_status.txt", 
            row.names = FALSE, 
            quote = FALSE, 
            sep = "\t")

# Print results
print(head(tau_results_combined))
print(summary_tau_combined)
print(wilcox_results_combined)

