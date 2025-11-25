### -Plot gene coverage between interacting gene pairs (updated, no CNEs)- ###
rm(list = ls())

library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(ggpubr)

# Define new file paths for interacting gene data
files_gene_pairs <- list(
  "eupsc_all_species_interacting" = "/Users/users/Desktop/Micro-C/topology_strength_analysis/eupsc_100k/409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof.eupsc_interacting_all_species_intervening_genes.txt",
  "eupsc_deca_only_interacting" = "/Users/users/Desktop/Micro-C/topology_strength_analysis/eupsc_100k/409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof.eupsc_interacting_deca_only_intervening_genes.txt",
  "eupsc_octbi_only_interacting" = "/Users/users/Desktop/Micro-C/topology_strength_analysis/eupsc_100k/409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof.eupsc_interacting_octbi_only_intervening_genes.txt",
  "eupsc_not_interacting" = "/Users/users/Desktop/Micro-C/topology_strength_analysis/eupsc_100k/409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof.eupsc_not_interacting_any_species_intervening_genes.txt",
  "octbi_all_species_interacting" = "/Users/users/Desktop/Micro-C/topology_strength_analysis/eupsc_100k/409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof.octbi_interacting_all_species_intervening_genes.txt",
  "octbi_deca_only_interacting" = "/Users/users/Desktop/Micro-C/topology_strength_analysis/eupsc_100k/409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof.octbi_interacting_deca_only_intervening_genes.txt",
  "octbi_octbi_only_interacting" = "/Users/users/Desktop/Micro-C/topology_strength_analysis/eupsc_100k/409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof.octbi_interacting_octbi_only_intervening_genes.txt",
  "octbi_not_interacting" = "/Users/users/Desktop/Micro-C/topology_strength_analysis/eupsc_100k/409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof.octbi_not_interacting_any_species_intervening_genes.txt"
)

# Load and format data
gene_data <- bind_rows(
  lapply(names(files_gene_pairs), function(name) {
    parts <- str_split(name, "_", simplify = TRUE)
    species <- parts[1]
    interaction_label <- paste(parts[2:ncol(parts)], collapse = "_")
    read_tsv(files_gene_pairs[[name]], show_col_types = FALSE) %>%
      mutate(
        species = species,
        interaction_detail = interaction_label,
        intervening_count = as.numeric(intervening_count),
        gene_coverage_percentage = as.numeric(gsub("%", "", gene_coverage_percentage))
      ) %>%
      filter(intervening_count > 0) %>% #At least one intervening gene must be present
      select(species, interaction_detail, gene_coverage_percentage)
  })
)

# Clean factors
gene_data$interaction_detail <- factor(gene_data$interaction_detail, levels = c(
  "all_species_interacting", "deca_only_interacting", "octbi_only_interacting", "not_interacting"
))

# Plot
p <- ggplot(gene_data, aes(x = interaction_detail, y = gene_coverage_percentage, fill = species)) +
  geom_boxplot(alpha = 0.8, outlier.size = 1, position = position_dodge(width = 0.8)) +
  theme_minimal() +
  labs(title = "Gene coverage percentage in gene pairs",
       x = "Interaction type",
       y = "Gene coverage (%)") +
  ylim(0, 100) +
  scale_fill_manual(values = c("#9FE2BF", "#6B5CAE"), name = "Species") +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1)) +
  scale_x_discrete(labels = c("All species", "Deca only", "Octbi only", "Not interacting"))

print(p)

ggsave("/Users/users/Desktop/Micro-C/figs_for_paper/Figure3/interacting_gene_coverage_box_4_cats.tiff", p, dpi = 200, width = 6, height = 6, units = "in")

# Wilcoxon tests
wilcox_results <- gene_data %>%
  group_by(species) %>%
  summarise(pw_test = list(pairwise.wilcox.test(gene_coverage_percentage, interaction_detail, p.adjust.method = "none")), .groups = "drop") %>%
  mutate(results = lapply(pw_test, function(pw) {
    res <- as.data.frame(as.table(pw$p.value))
    names(res) <- c("Group1", "Group2", "p_value")
    res <- res[!is.na(res$p_value), ]
    res$p_adjusted <- p.adjust(res$p_value, method = "BH")
    return(res)
  })) %>%
  select(species, results) %>%
  unnest(results)

print(wilcox_results)

write.table(wilcox_results, "/Users/users/Desktop/Micro-C/tables_for_paper/coverage/interacting_gene_coverage_wilcox_4cats.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Medians
gene_medians <- gene_data %>%
  group_by(species, interaction_detail) %>%
  summarise(median_gene_coverage = round(median(gene_coverage_percentage, na.rm = TRUE), 4))

print(gene_medians)

#COMPARE TO LOOPS
# Load loop files
loop_files <- list(
  eupsc = "/Users/users/Desktop/Micro-C/diff_loop_analysis/eupsc/output/eupsc_29cat_50k+100k/eupsc_29cat_50k+100k.tsv.intervening_genes",
  octbi = "/Users/users/Desktop/Micro-C/diff_loop_analysis/octbi/output/octbi_wt_50k+100k/octbi_wt_50k+100k.tsv.intervening_genes"
)

loop_data <- bind_rows(lapply(names(loop_files), function(sp) {
  read_tsv(loop_files[[sp]], show_col_types = FALSE) %>%
    mutate(
      species = sp,
      region_type = "loop",
      gene_coverage_percentage = as.numeric(gsub("%", "", gene_coverage_percentage))
    ) %>%
    filter(!is.na(gene_coverage_percentage), gene_coverage_percentage > 0) %>%
    select(species, region_type, gene_coverage_percentage)
}))

# Load existing gene pair data (from previous steps or file)
# Make sure this object exists already or replace with your actual data loading
# Example placeholder:
# gene_data <- read_tsv("your_cleaned_gene_pair_data.tsv")

# We'll assume 'gene_data' is already loaded with columns: species, interaction_detail, gene_coverage_percentage
# Add region_type column to distinguish gene pairs
gene_data$region_type <- paste0("gene_pair_", gene_data$interaction_detail)

# Combine datasets
combined_data <- bind_rows(
  loop_data,
  gene_data %>% select(species, region_type, gene_coverage_percentage)
)

# Plot comparison
combined_data$region_type <- factor(combined_data$region_type, levels = c(
  "gene_pair_all_species_interacting",
  "gene_pair_deca_only_interacting",
  "gene_pair_octbi_only_interacting",
  "gene_pair_not_interacting",
  "loop"
))

ggplot(combined_data, aes(x = region_type, y = gene_coverage_percentage, fill = species)) +
  geom_boxplot(alpha = 0.8) +
  theme_minimal() +
  labs(x = "Region type", y = "Gene coverage (%)", title = "Gene coverage in gene pairs and loops") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("eupsc" = "#9FE2BF", "octbi" = "#6B5CAE"))

# Wilcoxon tests: loop vs all gene pair categories (within each species)
categories_to_compare <- unique(combined_data$region_type)
categories_to_compare <- categories_to_compare[categories_to_compare != "loop"]

wilcox_loop_all_gp <- bind_rows(lapply(unique(combined_data$species), function(sp) {
  loop_values <- combined_data %>%
    filter(species == sp, region_type == "loop") %>%
    pull(gene_coverage_percentage)
  
  lapply(categories_to_compare, function(cat) {
    gp_values <- combined_data %>%
      filter(species == sp, region_type == cat) %>%
      pull(gene_coverage_percentage)
    
    result <- wilcox.test(loop_values, gp_values)
    
    tibble(
      species = sp,
      group1 = "loop",
      group2 = cat,
      p_value = result$p.value
    )
  }) %>% bind_rows()
})) %>%
  mutate(p_adjusted = p.adjust(p_value, method = "BH"))

# View results
print(wilcox_loop_all_gp)

# Save to file
write.table(wilcox_loop_all_gp, "/Users/users/Desktop/Micro-C/tables_for_paper/coverage/loop_vs_all_gene_pairs_wilcox.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Compute loop gene coverage medians
loop_medians <- loop_data %>%
  group_by(species) %>%
  summarise(
    region_type = "loop",
    median_gene_coverage = round(median(gene_coverage_percentage, na.rm = TRUE), 2),
    .groups = "drop"
  )

# View
print(loop_medians)




