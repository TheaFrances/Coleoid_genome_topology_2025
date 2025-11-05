library(dplyr)
library(ggplot2)
library(ggridges)

rm(list = ls())

# ---- Helper: load bedtools coverage output ----
# Expects 7 columns as described above
read_cov <- function(file) {
  df <- read.table(file, header = FALSE)
  stopifnot(ncol(df) >= 7)
  colnames(df)[1:7] <- c("chr","start","end","n_b_overlaps","bp_covered","a_len","frac_covered")
  df
}

# ---- Load per-stage coverage (anchor, non-anchor, interloop, non-interloop) ----
# Stage 20
anchor20      <- read_cov("cov_eupsc_20_anchor.tsv")
nonanchor20   <- read_cov("cov_eupsc_20_nonanchor.tsv")
interloop20   <- read_cov("cov_eupsc_20_interloop.tsv")
noninter20    <- read_cov("cov_eupsc_20_noninterloop.tsv")

# Stage 25
anchor25      <- read_cov("cov_eupsc_25_anchor.tsv")
nonanchor25   <- read_cov("cov_eupsc_25_nonanchor.tsv")
interloop25   <- read_cov("cov_eupsc_25_interloop.tsv")
noninter25    <- read_cov("cov_eupsc_25_noninterloop.tsv")

# Stage 29
anchor29      <- read_cov("cov_eupsc_29_anchor.tsv")
nonanchor29   <- read_cov("cov_eupsc_29_nonanchor.tsv")
interloop29   <- read_cov("cov_eupsc_29_interloop.tsv")
noninter29    <- read_cov("cov_eupsc_29_noninterloop.tsv")

# ---- Combine into long table using fraction covered ----
combined <- bind_rows(
  transform(anchor20,    value = frac_covered, group = "Anchor",        stage = "Stage 20"),
  transform(nonanchor20, value = frac_covered, group = "Non-Anchor",    stage = "Stage 20"),
  transform(interloop20, value = frac_covered, group = "Interloop",     stage = "Stage 20"),
  transform(noninter20,  value = frac_covered, group = "Non-Interloop", stage = "Stage 20"),
  
  transform(anchor25,    value = frac_covered, group = "Anchor",        stage = "Stage 25"),
  transform(nonanchor25, value = frac_covered, group = "Non-Anchor",    stage = "Stage 25"),
  transform(interloop25, value = frac_covered, group = "Interloop",     stage = "Stage 25"),
  transform(noninter25,  value = frac_covered, group = "Non-Interloop", stage = "Stage 25"),
  
  transform(anchor29,    value = frac_covered, group = "Anchor",        stage = "Stage 29"),
  transform(nonanchor29, value = frac_covered, group = "Non-Anchor",    stage = "Stage 29"),
  transform(interloop29, value = frac_covered, group = "Interloop",     stage = "Stage 29"),
  transform(noninter29,  value = frac_covered, group = "Non-Interloop", stage = "Stage 29")
) |>
  dplyr::select(chr, start, end, group, stage, value)

# ---- Stats (Wilcoxon on fractions) ----
pvals <- c(
  wilcox.test(subset(combined, stage=="Stage 20" & group=="Anchor")$value,
              subset(combined, stage=="Stage 20" & group=="Non-Anchor")$value)$p.value,
  wilcox.test(subset(combined, stage=="Stage 20" & group=="Interloop")$value,
              subset(combined, stage=="Stage 20" & group=="Non-Interloop")$value)$p.value,
  wilcox.test(subset(combined, stage=="Stage 20" & group=="Interloop")$value,
              subset(combined, stage=="Stage 20" & group=="Anchor")$value)$p.value,
  
  wilcox.test(subset(combined, stage=="Stage 25" & group=="Anchor")$value,
              subset(combined, stage=="Stage 25" & group=="Non-Anchor")$value)$p.value,
  wilcox.test(subset(combined, stage=="Stage 25" & group=="Interloop")$value,
              subset(combined, stage=="Stage 25" & group=="Non-Interloop")$value)$p.value,
  wilcox.test(subset(combined, stage=="Stage 25" & group=="Interloop")$value,
              subset(combined, stage=="Stage 25" & group=="Anchor")$value)$p.value,
  
  wilcox.test(subset(combined, stage=="Stage 29" & group=="Anchor")$value,
              subset(combined, stage=="Stage 29" & group=="Non-Anchor")$value)$p.value,
  wilcox.test(subset(combined, stage=="Stage 29" & group=="Interloop")$value,
              subset(combined, stage=="Stage 29" & group=="Non-Interloop")$value)$p.value,
  wilcox.test(subset(combined, stage=="Stage 29" & group=="Interloop")$value,
              subset(combined, stage=="Stage 29" & group=="Anchor")$value)$p.value
)
names(pvals) <- c(
  "20: Anchor vs Non-anchor", "20: Interloop vs Non-interloop", "20: Interloop vs Anchor",
  "25: Anchor vs Non-anchor", "25: Interloop vs Non-interloop", "25: Interloop vs Anchor",
  "29: Anchor vs Non-anchor", "29: Interloop vs Non-interloop", "29: Interloop vs Anchor"
)
adj <- p.adjust(pvals, method = "BH")
cat("Raw p-values:\n"); print(pvals)
cat("\nBH-adjusted p-values:\n"); print(adj)

# ---- Quick summaries ----
stats_summary <- combined |>
  group_by(stage, group) |>
  summarise(mean_frac = mean(value), median_frac = median(value), .groups = "drop")
print(stats_summary)

# ---- Plots ----
ggplot(combined, aes(x = value, y = group, fill = stage)) +
  geom_density_ridges(alpha = 0.6, scale = 1.2, rel_min_height = 0.01) +
  facet_wrap(~stage) +
  scale_fill_manual(values = c("Stage 20" = "goldenrod",
                               "Stage 25" = "#BB5566",
                               "Stage 29" = "#004488")) +
  xlab("Fraction of region covered by ATAC peaks") +
  ylab("") +
  theme_minimal(base_size = 13)

combined |>
  group_by(group, stage) |>
  summarise(mean = mean(value), se = sd(value)/sqrt(n()), .groups = "drop") |>
  ggplot(aes(x = group, y = mean, fill = stage)) +
  geom_col(position = position_dodge()) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(0.9), width = 0.2) +
  ylab("Mean fraction covered") +
  xlab("") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Fit a linear model (ANOVA-style) to test how ATAC coverage (value) depends on the type of 3D structure and stage
combined$group <- factor(combined$group,
                         levels = c("Non-Anchor","Anchor","Non-Interloop","Interloop"))
combined$stage <- factor(combined$stage, levels = c("Stage 20","Stage 25","Stage 29"))

lm_model <- lm(value ~ group * stage, data = combined)
summary(lm_model)
anova(lm_model)


