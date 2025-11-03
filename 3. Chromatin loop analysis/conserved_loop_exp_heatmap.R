# Plot conserved loop anchor gene expression - example for E. scolopes and S. officinalis conserved loops.

rm(list = ls())

# Load libraries
library(dplyr)
library(tidyr)
library(readr)
library(pheatmap)

# Load TPM matrix
expression_data <- read.delim("TPM_countDataMatrixNotFilt.txt", header = TRUE)
head(expression_data)

# Assign readable tissue names
tissue_names <- c("gene_id", "B1 (First arm left of hectocotylus)", "Mantle", "Left gill", "Hectocotylus (A1)", "Skin", 
                  "Testes", "Left optic lobe", "B4 Arm", "Right tentacle", "Right gill", "Supraesophogeal lobe", 
                  "Subesophogeal lobe", "Left optic lobe", "Right optic lobe", "Left white body", "Right white body", 
                  "Left gill", "Right gill", "Testes", "Subesophogeal lobe", "Left optic lobe", "Left white body", 
                  "Right white body", "Left gill", "Right gill", "Ovaries", "ANG", "Central core", "Central brain", 
                  "Left optic lobe", "Right optic lobe", "Left white body", "Right white body", "Ovaries", "ANG", 
                  "Central core", "ANG", "Central core", "ANG", "Central core", "ANG")

colnames(expression_data) <- make.names(tissue_names, unique = TRUE)

# Log-transform expression data
expression_data[,-1] <- lapply(expression_data[,-1], function(x) log(as.numeric(x) + 0.01))

# Average replicates per tissue
group_map <- list(
  "B1 arm" = c("B1..First.arm.left.of.hectocotylus."),
  "B4 arm" = c("B4.Arm"),
  "Tentacle" = c("Right.tentacle"),
  "Mantle" = c("Mantle"),
  "Hectocotylus" = c("Hectocotylus..A1."),
  "Skin" = c("Skin"),
  "Testes" = c("Testes", "Testes.1"),
  "Left optic lobe" = c("Left.optic.lobe", "Left.optic.lobe.1", "Left.optic.lobe.2", "Left.optic.lobe.3"),
  "Right optic lobe" = c("Right.optic.lobe", "Right.optic.lobe.1"),
  "Left gill" = c("Left.gill", "Left.gill.1", "Left.gill.2"),
  "Right gill" = c("Right.gill", "Right.gill.1", "Right.gill.2"),
  "Supraesophageal lobe" = c("Supraesophogeal.lobe"),
  "Subesophageal lobe" = c("Subesophogeal.lobe", "Subesophogeal.lobe.1"),
  "Left white body" = c("Left.white.body", "Left.white.body.1", "Left.white.body.2"),
  "Right white body" = c("Right.white.body", "Right.white.body.1", "Right.white.body.2"),
  "Ovaries" = c("Ovaries", "Ovaries.1"),
  "ANG" = c("ANG", "ANG.1", "ANG.2", "ANG.3", "ANG.4"),
  "Central core" = c("Central.core", "Central.core.1", "Central.core.2", "Central.core.3", "Central.core.4"),
  "Central brain" = c("Central.brain")
)

# Convert to long, average, then back to wide
averaged_data <- expression_data %>%
  pivot_longer(-gene_id, names_to = "sample", values_to = "logTPM") %>%
  mutate(tissue = sapply(sample, function(x) {
    match <- names(Filter(function(v) x %in% v, group_map))
    if (length(match) == 0) return(x) else return(match)
  })) %>%
  group_by(gene_id, tissue) %>%
  summarise(avg_logTPM = mean(logTPM, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = tissue, values_from = avg_logTPM)

# Function to extract cluster IDs from conserved loop files
extract_clusters <- function(filepath) {
  df <- read.delim(filepath, header = TRUE, stringsAsFactors = FALSE)
  
  # Extract just the two eupsc gene columns
  gene_cols <- df[, c("eupsc_genes_bin1", "eupsc_genes_bin2")]
  
  # Flatten, split on comma, clean up, and keep only cluster-prefixed values
  clusters <- gene_cols %>%
    unlist(use.names = FALSE) %>%
    as.character() %>%
    strsplit(",\\s*") %>%
    unlist() %>%
    trimws() %>%
    grep("^cluster_", ., value = TRUE) %>%
    unique()
  
  return(clusters)
}


# Load cluster IDs from conserved loop file
conserved_file <- "eupsc_sepof_consv_loops_50k+100k.txt"
head(conserved_file)

clusters_of_interest <- extract_clusters(conserved_file)

# Subset expression matrix to genes in conserved loops
heatmap_data <- averaged_data %>%
  filter(gene_id %in% clusters_of_interest) %>%
  column_to_rownames("gene_id")

# Plot heatmap with row scaling
p <- pheatmap(heatmap_data,
         scale = "row",
         cluster_rows = TRUE,
         treeheight_row = 0,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         fontsize_col = 10,
         border_color = NA,
         main = "Expression of E. scolopes and S. officinalis conserved loop genes")

# Save to file
ggsave(
  filename = "conserved_loop_expression_euspsc_sepof_tissues.tiff",
  plot = p,
  device = "tiff",
  width = 9,
  height = 6,
  units = "in",
  dpi = 200
)


