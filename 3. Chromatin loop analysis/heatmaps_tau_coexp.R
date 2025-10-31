#Diff loop heatmap clustered by LOOP in eupsc.
rm(list = ls())

# Load necessary libraries
library(dplyr)
library(pheatmap)
library(tidyverse)
library(readr)
library(stringr)
library(RColorBrewer)

#Explanation of this script----
#1. Merging Expression Data with Loop Mapping
#When you use inner_join(loop_gene_mapping, by = c("gene_id" = "genes")), genes that appear in multiple loops will be included in the merged_data for each loop they belong to.

#2. Sorting and Preparing Data for Heatmap
#After joining the data and sorting by loop_id, the sorted_data data frame will include all occurrences of genes across different loops. This means that genes appearing in multiple loops will be listed separately for each loop.

#3. Creating the Heatmap
#When converting sorted_expression_data to a matrix and then creating the heatmap with pheatmap, each row in sorted_expression_matrix corresponds to a unique gene-loop pair. If a gene is in multiple loops, it will have multiple rows in the matrix, one for each loop it appears in. The gaps_row parameter will create visual gaps based on the unique loop IDs, and rows with the same gene will be grouped correctly.

#4. Visualization
#The gaps_row parameter in pheatmap is used to create vertical lines (or gaps) between different loops. This helps in visually distinguishing the rows corresponding to different loops. Genes that are repeated across loops will be plotted separately for each loop they belong to.

# Function to read and process the loop gene mapping file----
read_loop_gene_mapping <- function(loop_file_path) {
  loop_data <- read_delim(loop_file_path, delim = "\t", col_names = TRUE)
  
  # Ensure that bin1_start and bin2_start are treated as full numbers
  loop_data <- loop_data %>%
    mutate(
      `bin1_start (bp)` = as.character(as.integer(`bin1_start (bp)`)),
      `bin2_start (bp)` = as.character(as.integer(`bin2_start (bp)`))
    ) %>%
    rowwise() %>%
    mutate(
      genes_bin1 = list(strsplit(genes_bin1, ",")[[1]]),
      genes_bin2 = list(strsplit(genes_bin2, ",")[[1]])
    ) %>%
    ungroup() %>%
    pivot_longer(cols = starts_with("genes_bin"), names_to = "bin", values_to = "gene_id") %>%
    unnest(gene_id) %>%
    mutate(gene_id = str_trim(gene_id)) %>%
    group_by(gene_id) %>%
    reframe(
      loop_id = paste(chromosome, min(`bin1_start (bp)`, `bin2_start (bp)`), max(`bin1_start (bp)`, `bin2_start (bp)`), sep = "_")
    )
  
  return(loop_data)
}

# Path to the loop gene mapping file
loop_file_path <- "tot_loops_eupsc_29_50k+100k.tsv.genes_rm_dups"

# Read the loop gene mapping
loop_gene_mapping <- read_loop_gene_mapping(loop_file_path)

# Remove duplicate lines (same gene and same loop)
loop_gene_mapping <- loop_gene_mapping %>%
  distinct()

#Print original number of loops
cat("Original loops (stage 29):", n_distinct(loop_gene_mapping$loop_id), "\n")

# Read expression data----
genes_to_search_29 <- readLines("tot_loops_eupsc_29_50k+100k.tsv.genes_list.txt")

tpm_data <- read.table("/Users/users/Desktop/Micro-C/expression_data/Hannah_TPM_normalizedcounts_development.txt", 
                       header = TRUE, sep = "\t")

# Filter the data for each stage
pattern29 <- paste0("\\b(", paste(genes_to_search_29, collapse = "|"), ")\\b")
expression_data <- tpm_data[grep(pattern29, tpm_data$gene_id), ]

# Replace NA values with 0 (if there are any)
expression_data[,-1][is.na(expression_data[,-1])] <- 0

# Remove rows that have 0 expression for every tissue
expression_data_0_rm <- expression_data[rowSums(expression_data[,-1] == 0) != (ncol(expression_data) - 1), ]

# Log-transform data
expression_data_0_rm[,-1] <- lapply(expression_data_0_rm[,-1], function(x) {
  log(as.numeric(as.character(x)) + 0.01)
})

# Load loop-gene mappings and tag stage
loop_20 <- read_loop_gene_mapping("tot_loops_eupsc_20_50k+100k.tsv.genes_rm_dups") %>% mutate(stage = "20")
loop_25 <- read_loop_gene_mapping("tot_loops_eupsc_25_50k+100k.tsv.genes_rm_dups") %>% mutate(stage = "25")
loop_29 <- loop_gene_mapping %>% mutate(stage = "29")

# Combine and summarize gene sets
all_loops <- bind_rows(loop_20, loop_25, loop_29)

loop_signatures <- all_loops %>%
  group_by(stage, loop_id) %>%
  summarise(gene_set = paste(sort(unique(gene_id)), collapse = ","), .groups = "drop")

# Identify gene sets used in multiple stages
shared_sets <- loop_signatures %>%
  add_count(gene_set) %>%
  group_by(gene_set) %>%
  filter(n_distinct(stage) > 1) %>%
  pull(gene_set) %>% 
  unique()

# Keep only stage 29 loops with unique gene sets
loop_ids_to_keep <- loop_signatures %>%
  filter(stage == "29", !(gene_set %in% shared_sets)) %>%
  pull(loop_id)

# Filter loop-gene mapping
loop_gene_mapping <- loop_gene_mapping %>%
  filter(loop_id %in% loop_ids_to_keep)

# Merge expression_data with loop_gene_mapping
merged_data <- expression_data_0_rm %>%
  inner_join(loop_gene_mapping, by = "gene_id")

# Sort merged data by loop_id
# Ensure loop_id is sorted numerically based on the numeric portion of the ID
sorted_data <- merged_data %>%
  arrange(as.numeric(str_extract(loop_id, "\\d+")))

# Filter out loops with only one gene
loop_gene_count <- sorted_data %>%
  group_by(loop_id) %>%
  summarise(gene_count = n(), .groups = 'drop')

loops_with_multiple_genes <- loop_gene_count %>%
  filter(gene_count > 1) %>%
  pull(loop_id)

filtered_sorted_data <- sorted_data %>%
  filter(loop_id %in% loops_with_multiple_genes)

# Update Loop_ID names in row_annotation
row_annotation <- data.frame(Loop_ID = filtered_sorted_data$loop_id)
rownames(row_annotation) <- filtered_sorted_data$gene_id

# Replace Loop_ID prefixes dynamically and remove unwanted suffixes
row_annotation$Loop_ID <- sapply(
  row_annotation$Loop_ID,
  function(x) {
    group_number <- as.numeric(gsub(".*Lachesis_group(\\d+).*", "\\1", x)) + 1
    gsub("Lachesis_group\\d+__.*?length_\\d+_", paste0("Chromosome_", group_number, "_"), x)
  }
)

# Ensure expression data is numeric
sorted_expression_data <- filtered_sorted_data %>%
  dplyr::select(-loop_id) %>%
  dplyr::select(-gene_id) %>%
  mutate(across(everything(), as.numeric))

# Convert to matrix for pheatmap
sorted_expression_matrix <- as.matrix(sorted_expression_data)
rownames(sorted_expression_matrix) <- filtered_sorted_data$gene_id

# Recalculate row gaps based on sorted Loop_IDs
rle_results <- rle(as.character(row_annotation$Loop_ID))
row_gaps_filtered <- cumsum(rle_results$lengths)

# Define colors for Loop_IDs
large_color_palette <- colorRampPalette(c(
  brewer.pal(8, "Set3"),
  brewer.pal(4, "Dark2"),
  brewer.pal(4, "Paired"),
  brewer.pal(7, "Spectral"),
  brewer.pal(5, "Set1")
))(length(unique(row_annotation$Loop_ID)))

annotation_colors <- list(
  Loop_ID = setNames(large_color_palette, unique(row_annotation$Loop_ID))
)

sorted_expression_matrix_3stages <- sorted_expression_matrix[,!grepl("^st14_|^st16_|^st18_|^st22_", colnames(sorted_expression_matrix))
]


# Generate the heatmap
heatmap_loops_1genes_rm_with_labels_3stages <- pheatmap(
  sorted_expression_matrix_3stages, 
  cluster_rows = FALSE,  
  cluster_cols = FALSE, 
  scale = 'row', 
  show_rownames = FALSE, 
  border_color = NA,
  gaps_row = row_gaps_filtered,
  cellwidth = 10,        
  cellheight = 3.5,      
  fontsize = 8,          
  width = 2.5,           
  height = 34,           
  annotation_row = row_annotation, 
  annotation_colors = annotation_colors
)

# Save the heatmap to files
ggsave(heatmap_loops_1genes_rm_with_labels_3stages, 
       file = "heatmap_stage_29_only_loop_genes_50k_100k.tiff", 
       height = 9.5, width = 7, limitsize = FALSE)

#Print number of loops after filtering
cat("Loops kept after gene set filtering:", length(loop_ids_to_keep), "\n")


#Tau
calculate_tau <- function(expression_values) {
  x_max <- max(expression_values, na.rm = TRUE)
  if (x_max <= 0) return(0)  # Avoid division by zero or negative max
  n <- sum(!is.na(expression_values))  # Count only non-NA values
  if (n <= 1) return(0)  # If only one tissue, Tau is 0
  
  tau <- sum(1 - (expression_values / x_max), na.rm = TRUE) / (n - 1)
  return(max(0, min(tau, 1)))  # Ensure Tau is between 0 and 1
}
# Tau on log-transformed expression data

# Use log-transformed expression matrix
log_expr_data <- expression_data_0_rm  # already log(TPM + 0.01)

# Filter to loop-relevant genes
log_expr_data <- log_expr_data %>%
  filter(gene_id %in% loop_gene_mapping$gene_id)

# Join with loop info
expression_with_loop <- inner_join(log_expr_data, loop_gene_mapping, by = "gene_id")

# Compute Tau per gene
expression_with_loop <- expression_with_loop %>%
  rowwise() %>%
  mutate(tau = calculate_tau(c_across(where(is.numeric)))) %>%
  ungroup()

# Summarize Tau by loop
tau_by_loop <- expression_with_loop %>%
  group_by(loop_id) %>%
  summarise(mean_tau = mean(tau, na.rm = TRUE),
            median_tau = median(tau, na.rm = TRUE),
            gene_count = n())

# Overall Tau summary
overall_mean_tau <- mean(tau_by_loop$mean_tau, na.rm = TRUE)
overall_median_tau <- median(tau_by_loop$median_tau, na.rm = TRUE)

cat("Overall mean of mean_tau across loops (log-transformed):", overall_mean_tau, "\n")
cat("Overall median of median_tau across loops (log-transformed):", overall_median_tau, "\n")


library(tidyr)

# Function to calculate mean pairwise correlation within a loop
calculate_loop_coexpression <- function(expr_matrix) {
  if (nrow(expr_matrix) < 2) return(NA)  # Cannot correlate fewer than 2 genes
  cor_matrix <- cor(t(expr_matrix), use = "pairwise.complete.obs")
  upper_tri_values <- cor_matrix[upper.tri(cor_matrix)]
  return(mean(upper_tri_values, na.rm = TRUE))
}

# Prepare expression matrix: rows = genes, columns = tissues
loop_expr_long <- log_expr_data %>%
  pivot_longer(-gene_id, names_to = "tissue", values_to = "expression")

# Add loop info
loop_expr_long <- inner_join(loop_expr_long, loop_gene_mapping, by = "gene_id")

# Spread to wide format per loop
coexpression_by_loop <- loop_expr_long %>%
  pivot_wider(id_cols = c(loop_id, gene_id), names_from = tissue, values_from = expression) %>%
  group_by(loop_id) %>%
  group_modify(~ {
    expr_matrix <- as.data.frame(.x[,-(1:2)])  # remove loop_id and gene_id
    coexp <- calculate_loop_coexpression(expr_matrix)
    tibble(coexpression = coexp)
  }) %>%
  ungroup()

# Merge with Tau
tau_coexp_by_loop <- tau_by_loop %>%
  left_join(coexpression_by_loop, by = "loop_id")

# Save or print
cat("Mean coexpression across loops:", mean(tau_coexp_by_loop$coexpression, na.rm = TRUE), "\n")
cat("Median coexpression across loops:", median(tau_coexp_by_loop$coexpression, na.rm = TRUE), "\n")


