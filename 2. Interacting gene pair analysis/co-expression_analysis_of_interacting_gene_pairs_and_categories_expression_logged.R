# Co-expression analysis of gene pairs in different interaction categories in E. scolopes, the species with available expression data.
# This code calculates Pearson’s correlation coefficients for co-expression per gene pair across E. scolopes tissues for each interaction category and plots it as a density plot.
# It also calculates significant differences (Wilcoxon test), means, and medians in co-expression for gene pairs across different interaction categories.

# Clear workspace
rm(list = ls())

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(broom)

# Reading the (TPM normalised) gene expression file of expression across E. scolopes tissues
expression_data <- read.delim("TPM_countDataMatrixNotFilt.txt", header=T)

# Define the new tissue names in the same order as the sample IDs
tissue_names <- c("gene_id", "B1 (First arm left of hectocotylus)", "Mantle", "Left gill", "Hectocotylus (A1)", "Skin", 
                  "Testes", "Left optic lobe", "B4 Arm", "Right tentacle", "Right gill", "Supraesophogeal lobe", 
                  "Subesophogeal lobe", "Left optic lobe", "Right optic lobe", "Left white body", "Right white body", 
                  "Left gill", "Right gill", "Testes", "Subesophogeal lobe", "Left optic lobe", "Left white body", 
                  "Right white body", "Left gill", "Right gill", "Ovaries", "ANG", "Central core (CC)", "Central Brain", 
                  "Left optic lobe", "Right optic lobe", "Left white body", "Right white body", "Ovaries", "ANG", 
                  "Central core", "ANG", "Central core", "ANG", "Central core", "ANG")

# Assign the column names to the data frame
colnames(expression_data) <- tissue_names

# Make column names unique otherwise you will have problems downstream
make_unique <- function(names) {
  unique_names <- make.names(names, unique = TRUE)
  return(unique_names)
}
colnames(expression_data) <- make_unique(colnames(expression_data))

# Verify the column names are unique
print(colnames(expression_data))

# Convert all columns except 'gene_id' to numeric and log-transform (adding a small pseudocount to avoid log(0))
expression_data[,-1] <- lapply(expression_data[,-1], function(x) {
  log(as.numeric(as.character(x)) + 0.01)
})

# Load interaction data (note used to be interactions_pecten_chr_genom_dist_ave_threshold10.txt which I think is the same file but renamed)
interaction_data <- read.delim("409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof.txt")
head(interaction_data)

# Remove duplicate rows (if there are any)
interaction_data<- distinct(interaction_data)

# Split the gene pairs into separate columns
interaction_data <- interaction_data %>%
  separate(Orth_pair_based_on_eupsc_interaction_matrix, into = c("Gene1", "Gene2"), sep = ",")

# Further split to isolate IDs for the species you have expression data for
interaction_data <- interaction_data %>%
  mutate(Gene1 = sub(";.*", "", Gene1),
         Gene2 = sub(";.*", "", Gene2))

head(interaction_data)

# Merge gene expression data with interaction data
merged_data <- interaction_data %>%
  inner_join(expression_data, by = c("Gene1" = "gene_id")) %>%
  inner_join(expression_data, by = c("Gene2" = "gene_id"), suffix = c(".Gene1", ".Gene2"))

head(merged_data)

# Function to calculate correlation between two genes across tissues
calculate_correlation <- function(row) {
  gene1_expression <- as.numeric(row[grep(".Gene1$", names(row))])
  gene2_expression <- as.numeric(row[grep(".Gene2$", names(row))])
  
  if (sd(gene1_expression) == 0 || sd(gene2_expression) == 0) {
    return(NA)  # Return NA if the standard deviation is zero
  } else {
    return(cor(gene1_expression, gene2_expression, use = "complete.obs"))
  }
}

#Note there are non finite values when there is zero variance: 
#If the expression values of one or both genes are constant across all tissues the standard deviation is zero. 
#The correlation is not defined when the standard deviation of one or both variables is zero, resulting in NA values which will not be plotted.
# NA values are also removed from the mean, median and Wilcox calculations.

# Apply the correlation calculation
merged_data <- merged_data %>%
  rowwise() %>%
  mutate(Correlation = calculate_correlation(cur_data()))

# Define the custom colors for filling
interaction_cols <- c(
  "interacting_all_species" = adjustcolor("#FF7878", alpha.f = 0.7), 
  "interacting_deca_only" = adjustcolor("palegreen", alpha.f = 0.7), 
  "interacting_octbi_only" = adjustcolor("mediumpurple1", alpha.f = 0.7), 
  "not_interacting_any_species" = adjustcolor("orange", alpha.f = 0.7)
)

# Plot correlation distributions
coexp_density <- ggplot(merged_data) +
  geom_density(aes(x = Correlation, fill = Interaction_status), alpha = 0.7) +
  labs(title = "Correlation of gene expression in different interaction categories",
       x = expression(paste("Correlation of gene expression (logTPM) in ", italic("E. scolopes"))),
       y = "Density") +
  scale_fill_manual(name = "Interaction status",
                    values = c(interaction_cols)) + theme_minimal() + xlim(-1,1) + ylim(0, 1.7) +
  theme(
    plot.title = element_text(size = 20),      # Title text size
    axis.title.x = element_text(size = 20),    # X-axis title text size
    axis.title.y = element_text(size = 20),    # Y-axis title text size
    axis.text.x = element_text(size = 20),     # X-axis label text size
    axis.text.y = element_text(size = 20)      # Y-axis label text size
  )

# Print the plot
print(coexp_density)

#NOTE 
#ggplot2 automatically extends the axis range slightly beyond your data to avoid cutting off density curves.
#The density function estimates a smooth distribution and might extend past the actual data range.

# Save the plot
ggsave(coexp_density, file = "coexp_density_all_cats_eupsc_with_sof_logged.tiff", height = 6, width = 10)

#STATS for density plots----
# Get unique interaction statuses
interaction_statuses <- unique(merged_data$Interaction_status)

# Calculate means and medians for each group
summary_stats <- merged_data %>%
  group_by(Interaction_status) %>%
  summarise(
    mean_correlation = mean(Correlation, na.rm = TRUE),
    median_correlation = median(Correlation, na.rm = TRUE)
  )

# Merge summary statistics back into merged_data
merged_data <- merged_data %>%
  left_join(summary_stats, by = "Interaction_status")

# Generate pairwise combinations of interaction statuses
pairwise_combinations <- combn(interaction_statuses, 2, simplify = FALSE)

# Function to perform Wilcoxon test
perform_wilcox_test <- function(pair) {
  group1_mean <- summary_stats %>% filter(Interaction_status == pair[1]) %>% pull(mean_correlation)
  group1_median <- summary_stats %>% filter(Interaction_status == pair[1]) %>% pull(median_correlation)
  group2_mean <- summary_stats %>% filter(Interaction_status == pair[2]) %>% pull(mean_correlation)
  group2_median <- summary_stats %>% filter(Interaction_status == pair[2]) %>% pull(median_correlation)
  
  group1 <- merged_data %>% filter(Interaction_status == pair[1]) %>% pull(Correlation)
  group2 <- merged_data %>% filter(Interaction_status == pair[2]) %>% pull(Correlation)
  
  test_result <- wilcox.test(group1, group2, paired = FALSE)
  
  tidy(test_result) %>%
    mutate(
      group1 = pair[1], 
      group2 = pair[2],
      mean_group1 = group1_mean,
      median_group1 = group1_median,
      mean_group2 = group2_mean,
      median_group2 = group2_median
    )
}

# Perform Wilcoxon tests for all pairwise combinations
wilcox_results <- pairwise_combinations %>%
  lapply(perform_wilcox_test) %>%
  bind_rows()

# Apply Benjamini-Hochberg correction for multiple testing
wilcox_results <- wilcox_results %>%
  mutate(adj_p_value = p.adjust(p.value, method = "BH"))  # Adjust p-values using BH correction

# Print the results in a table and save
print(wilcox_results)

write.table(wilcox_results, "coexp_eupsc_wilcox_with_sof_logged.txt", append = FALSE, sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE, quote = FALSE)

