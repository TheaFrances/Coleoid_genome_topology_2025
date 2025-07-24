# GO analyses for genes in gene pairs interacting across the coleoids category

# Clear workspace
rm(list = ls())

# Load required libraries
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(clusterProfiler)
library(org.Ooctbi.eg.db)
library(GO.db)

# Read interaction data
eupsc_100k_obi_50k_pecten_chr <- read.delim("409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof.txt")

# Subset to gene pairs interacting across coleoids
eupsc_100k_obi_50k_pecten_chr <- subset(eupsc_100k_obi_50k_pecten_chr, Interaction_status == "interacting_all_species")

# Extract unique O. bimaculoides gene IDs from ortholog pairs
interacting_all_species_octbi_genes <- eupsc_100k_obi_50k_pecten_chr %>%
  distinct(Orth_pair_based_on_eupsc_interaction_matrix) %>%
  mutate(octbi_id = str_extract_all(Orth_pair_based_on_eupsc_interaction_matrix, "XP_\\d+\\.\\d+")) %>%
  pull(octbi_id) %>%
  unlist()

# Save to file (optional)
writeLines(interacting_all_species_octbi_genes, "interacting_all_spp_octbi_ids.txt")

# Prepare gene list and background
file <- "interacting_all_spp_octbi_ids.txt"
genes <- read.table(file)
gs <- unique(as.character(genes[,1]))

# Background gene list
allgenes <- read.table("term2gene")
ag <- unique(as.character(allgenes[,2]))

# GO term names
go2name <- select(GO.db, keys = keys(GO.db), columns = columns(GO.db)[1:2])

# Run GO enrichment and plot
for (ont in c("CC", "MF", "BP")) {
  message("Running enrichment for: ", ont)
  
  ego <- enrichGO(
    gene          = gs,
    universe      = ag,
    keyType       = "GENENAME",
    OrgDb         = org.Ooctbi.eg.db,
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  write.table(ego, paste0(file, "_enrich.", ont, ".go"))
  
  if (nrow(ego) > 0) {
    p <- dotplot(ego, showCategory = 10) +
      theme(
        text = element_text(size = 8),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        plot.title = element_text(size = 10)
      ) +
      guides(size = guide_legend(override.aes = list(size = 3)))
    
    ggsave(paste0(file, "_enrich_dotplot_", ont, ".pdf"), plot = p, width = 6, height = 6, dpi = 200)
  }
}
