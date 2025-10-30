# Clear workspace
rm(list = ls())

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsignif)  # For significance annotations
library(gridExtra)  # For combining plots

# Clear the workspace
rm(list = ls())

# Read the list of genes for each stage
genes_to_search_20 <- readLines("/Users/users/Desktop/Micro-C/diff_loop_analysis/eupsc/output/20_tot_loop_count/tot_loops_20_50k+100k_method1.tsv.genes_rm_dups_list.txt")
genes_to_search_25 <- readLines("/Users/users/Desktop/Micro-C/diff_loop_analysis/eupsc/output/25_tot_loop_count/tot_loops_25_50k+100k_method1.tsv.genes_rm_dups_list.txt")
genes_to_search_29 <- readLines("/Users/users/Desktop/Micro-C/diff_loop_analysis/eupsc/output/29_tot_loop_count/tot_loops_29_50k+100k_method1.tsv.genes_rm_dups_list.txt")
#genes_to_search_all <- readLines("/Users/users/Desktop/Micro-C/diff_loop_analysis/eupsc/output/loops_in_all_stages_50k+100k.genes_list.txt")
#genes_to_search_29cat <- readLines("/Users/users/Desktop/Micro-C/diff_loop_analysis/eupsc/output/stage_29cat_only_loops.genes_list.txt")

# Read the TPM normalized data
tpm_data <- read.table("/Users/users/Desktop/Micro-C/expression_data/Hannah_TPM_normalizedcounts_development.txt", 
                       header = TRUE, sep = "\t")

# Average replicates for each developmental stage
stage_means <- tpm_data %>%
  rowwise() %>%
  mutate(
    Stage_14 = mean(c_across(starts_with("st14_")), na.rm = TRUE),
    Stage_16 = mean(c_across(starts_with("st16_")), na.rm = TRUE),
    Stage_18 = mean(c_across(starts_with("st18_")), na.rm = TRUE),
    Stage_20 = mean(c_across(starts_with("st20_")), na.rm = TRUE),
    Stage_22 = mean(c_across(starts_with("st22_")), na.rm = TRUE),
    Stage_24 = mean(c_across(starts_with("st24_")), na.rm = TRUE),
    Stage_27 = mean(c_across(starts_with("st27_")), na.rm = TRUE),
    Stage_28 = mean(c_across(starts_with("st28_")), na.rm = TRUE)
  ) %>%
  dplyr::select(gene_id, Stage_14, Stage_16, Stage_18, Stage_20, Stage_22, Stage_24, Stage_27, Stage_28)

# Filter the data for each stage
pattern20 <- paste0("\\b(", paste(genes_to_search_20, collapse = "|"), ")\\b")
TPM_norm_counts20 <- stage_means[grep(pattern20, stage_means$gene_id), ]

pattern25 <- paste0("\\b(", paste(genes_to_search_25, collapse = "|"), ")\\b")
TPM_norm_counts25 <- stage_means[grep(pattern25, stage_means$gene_id), ]

pattern29 <- paste0("\\b(", paste(genes_to_search_29, collapse = "|"), ")\\b")
TPM_norm_counts29 <- stage_means[grep(pattern29, stage_means$gene_id), ]

#patternall <- paste0("\\b(", paste(genes_to_search_all, collapse = "|"), ")\\b")
#TPM_norm_countsall <- stage_means[grep(patternall, stage_means$gene_id), ]

#pattern29cat <- paste0("\\b(", paste(genes_to_search_29cat, collapse = "|"), ")\\b")
#TPM_norm_counts29cat <- stage_means[grep(pattern29cat, stage_means$gene_id), ]

# Combine all datasets for plotting
TPM_norm_counts20$Stage <- "Stage 20"
TPM_norm_counts25$Stage <- "Stage 25"
TPM_norm_counts29$Stage <- "Stage 29"
#TPM_norm_countsall$Stage <- "All stages"
#TPM_norm_counts29cat$Stage <- "Stage 29cat"

combined_data <- bind_rows(TPM_norm_counts20, TPM_norm_counts25, TPM_norm_counts29) #, TPM_norm_countsall) #Add TPM_norm_counts29cat, TPM_norm_counts29all if you want to check that as well

#Remove genes that appear more than once, this should remove genes that appear in multiple stages
filtered_data <- combined_data %>%
  group_by(gene_id) %>%
  filter(n() == 1) %>%
  ungroup()

# Or not
#filtered_data <- combined_data

# Continue with your reshaping and plotting as before
filtered_data_long <- filtered_data %>%
  pivot_longer(cols = -c(gene_id, Stage), names_to = "Tissue", values_to = "Expression")

# Log-transform the expression data with a pseudocount of 0.01
filtered_data_long <- filtered_data_long %>%
  mutate(Log_Expression = log(Expression + 0.01))


# Generate all pairwise comparisons for significance testing
pairwise_comparisons <- combn(unique(filtered_data_long$Stage), 2, simplify = FALSE)

custom_colors <- c("Stage 20" = "goldenrod", "Stage 25" = "#BB5566", "Stage 29" = "#004488",  "All stages" = "grey") #Add 29cat if you want to check that as well

# Plot the boxplots with significance testing
boxplot <- ggplot(filtered_data_long, aes(x = Stage, y = Log_Expression, fill = Stage)) +
  geom_boxplot(outlier.color = "black", outlier.size = 1, alpha = 0.8) +
  facet_wrap(~ Tissue, scales = "free", ncol = 4) +  # Adjust to 4 plots per row
  labs(title = "Gene expression across development for developmental stage-specific loops",
       x = "Developmental stage",
       y = "Expression (logTPM)") +
  scale_fill_manual(values = custom_colors) +  # Use scale_fill_manual for custom fill colors
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_signif(
    comparisons = pairwise_comparisons,
    step_increase = 0.1,
    test = "wilcox.test",  # Use Wilcoxon test
    map_signif_level = TRUE
  )

# Display the plot
print(boxplot)

#Save
ggsave(boxplot, file = "/Users/users/Desktop/Micro-C/figs_for_paper/Supplementary/diff_loop_figs/boxplots_eupsc_dev_stage_specific_gene_exp_50k+100k.tiff", height = 9, width = 13)


#Extra----
#Print exact signigicance values:
# Create a data frame to store the significance test results
significance_results <- data.frame(
  Tissue = character(),
  Comparison = character(),
  Statistic = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each tissue to perform pairwise tests
unique_tissues <- unique(filtered_data_long$Tissue)
for (tissue in unique_tissues) {
  tissue_data <- filtered_data_long %>%
    filter(Tissue == tissue)
  
  # Perform pairwise comparisons for the stages
  comparisons <- combn(unique(tissue_data$Stage), 2, simplify = FALSE)
  for (comparison in comparisons) {
    stage1 <- comparison[1]
    stage2 <- comparison[2]
    
    # Perform Wilcoxon test
    test_result <- wilcox.test(
      tissue_data$Log_Expression[tissue_data$Stage == stage1],
      tissue_data$Log_Expression[tissue_data$Stage == stage2],
      exact = FALSE
    )
    
    # Add results to the data frame
    significance_results <- rbind(
      significance_results,
      data.frame(
        Tissue = tissue,
        Comparison = paste(stage1, "vs", stage2),
        Statistic = test_result$statistic,
        P_value = test_result$p.value
      )
    )
  }
}


# Combine all datasets for plotting
TPM_norm_counts20$Stage <- "Stage 20"
TPM_norm_counts25$Stage <- "Stage 25"
TPM_norm_counts29$Stage <- "Stage 29"
TPM_norm_counts29cat$Stage <- "Stage 29cat"

combined_data <- bind_rows(TPM_norm_counts20, TPM_norm_counts25, TPM_norm_counts29) #Add TPM_norm_counts29cat if you want to check that as well


# Identify genes that appear in multiple stages
genes_in_multiple_stages <- combined_data %>%
  group_by(gene_id) %>%
  summarise(stage_count = n_distinct(Stage)) %>%
  filter(stage_count > 1) %>%
  pull(gene_id)

# Filter out genes that appear in multiple stages
filtered_data <- combined_data %>%
  filter(!gene_id %in% genes_in_multiple_stages)

# Continue with your reshaping and plotting as before
filtered_data_long <- filtered_data %>%
  pivot_longer(cols = -c(gene_id, Stage), names_to = "Tissue", values_to = "Expression") %>%
  rename(Gene = gene_id)

# Log-transform the expression data with a pseudocount of 0.01
filtered_data_long <- filtered_data_long %>%
  mutate(Log_Expression = log(Expression + 0.01))


# Generate all pairwise comparisons for significance testing
pairwise_comparisons <- combn(unique(filtered_data_long$Stage), 2, simplify = FALSE)

custom_colors <- c("Stage 20" = "goldenrod", "Stage 25" = "#BB5566", "Stage 29" = "#004488") #Add 29cat if you want to check that as well


# Plot the boxplots with significance testing
boxplot <- ggplot(filtered_data_long, aes(x = Stage, y = Log_Expression, fill = Stage)) +
  geom_boxplot(outlier.color = "black", outlier.size = 1, alpha = 0.8) +
  facet_wrap(~ Tissue, scales = "free") +
  labs(title = "Gene expression across development for developmental-stage specific loops",
       x = "Developmental stage",
       y = "Expression (logTPM)") +
  scale_fill_manual(values = custom_colors) +  # Use scale_fill_manual for custom fill colors
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_signif(
    comparisons = pairwise_comparisons,
    step_increase = 0.1,
    test = "wilcox.test",  # Use Wilcoxon test
    map_signif_level = TRUE
  )

# Display the plot
print(boxplot)


d#Extra----
#Print exact significance values and means and medians
# Create a data frame to store the significance test results
significance_results <- data.frame(
  Tissue = character(),
  Comparison = character(),
  Statistic = numeric(),
  P_value = numeric(),
  Mean_Stage1 = numeric(),
  Median_Stage1 = numeric(),
  Mean_Stage2 = numeric(),
  Median_Stage2 = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each tissue to perform pairwise tests
unique_tissues <- unique(filtered_data_long$Tissue)
for (tissue in unique_tissues) {
  tissue_data <- filtered_data_long %>%
    filter(Tissue == tissue)
  
  # Perform pairwise comparisons for the stages
  comparisons <- combn(unique(tissue_data$Stage), 2, simplify = FALSE)
  for (comparison in comparisons) {
    stage1 <- comparison[1]
    stage2 <- comparison[2]
    
    # Subset data for each stage
    stage1_data <- tissue_data$Log_Expression[tissue_data$Stage == stage1]
    stage2_data <- tissue_data$Log_Expression[tissue_data$Stage == stage2]
    
    # Perform Wilcoxon test
    test_result <- wilcox.test(
      stage1_data,
      stage2_data,
      exact = FALSE
    )
    
    # Calculate mean and median for each stage
    mean_stage1 <- mean(stage1_data, na.rm = TRUE)
    median_stage1 <- median(stage1_data, na.rm = TRUE)
    mean_stage2 <- mean(stage2_data, na.rm = TRUE)
    median_stage2 <- median(stage2_data, na.rm = TRUE)
    
    # Add results to the data frame
    significance_results <- rbind(
      significance_results,
      data.frame(
        Tissue = tissue,
        Comparison = paste(stage1, "vs", stage2),
        Statistic = test_result$statistic,
        P_value = test_result$p.value,
        Mean_Stage1 = mean_stage1,
        Median_Stage1 = median_stage1,
        Mean_Stage2 = mean_stage2,
        Median_Stage2 = median_stage2
      )
    )
  }
}

# View or save the results
print(significance_results)
#write.csv(significance_results, "significance_results_with_means_medians.csv", row.names = FALSE)
