rm(list = ls())

# Load overlap count data
# This section loads loop, interloop, and background regions for E. scolopes,
# calculates region lengths, and normalizes CNE counts per kilobase.

library(dplyr)
library(ggplot2)

# Helper function: read overlap count files and normalise by region length (per kb)
read_overlap <- function(file) {
  df <- read.table(file, header = FALSE)
  colnames(df)[1:4] <- c("chr", "start", "end", "cne_overlaps")
  df$length <- df$end - df$start
  df$per_kb <- df$cne_overlaps / (df$length / 1000)
  return(df)
}
# Load data
eupsc_anchor <- read_overlap("eupsc_loop_anchors_CNE_overlap_counts.bed")
eupsc_nonanchor <- read_overlap("eupsc_non_anchor_CNE_overlap_counts.bed")
eupsc_interloop <- read_overlap("eupsc_interloop_CNE_overlap_counts.bed")
eupsc_noninterloop <- read_overlap("eupsc_non_interloop_CNE_overlap_counts.bed")

# Combine into a single dataframe for plotting
eupsc_anchor$type <- "Anchor"
eupsc_nonanchor$type <- "Non-Anchor"
eupsc_interloop$type <- "Interloop"
eupsc_noninterloop$type <- "Non-Interloop"

combined_eupsc <- bind_rows(
  eupsc_anchor,
  eupsc_nonanchor,
  eupsc_interloop,
  eupsc_noninterloop
)

# Ensure correct factor order for plotting
combined_eupsc$type <- factor(combined_eupsc$type,
                              levels = c("Non-Anchor", "Anchor", "Non-Interloop", "Interloop"))

# Plot: Boxplot of CNE overlap per kb
ggplot(combined_eupsc, aes(x = type, y = per_kb, fill = type)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  xlab("") +
  ylab("CNE overlaps per kb") +
  scale_fill_brewer(palette = "Set2")

# Statistical comparisons (Wilcoxon rank-sum tests)
pvals <- c(
  wilcox.test(eupsc_anchor$per_kb, eupsc_nonanchor$per_kb)$p.value,
  wilcox.test(eupsc_interloop$per_kb, eupsc_noninterloop$per_kb)$p.value,
  wilcox.test(eupsc_interloop$per_kb, eupsc_anchor$per_kb)$p.value
)

names(pvals) <- c(
  "Anchor vs Non-anchor",
  "Interloop vs Non-interloop",
  "Interloop vs Anchor"
)

# Adjust with BH correction
pvals_adj <- p.adjust(pvals, method = "BH")

# Print results
cat("Raw p-values:\n")
print(round(pvals, 15))

cat("\nBH-adjusted p-values:\n")
print(round(pvals_adj, 15))
