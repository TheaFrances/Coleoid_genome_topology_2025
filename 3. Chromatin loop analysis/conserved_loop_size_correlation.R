rm(list = ls())

# Load libraries
library(tidyr)
library(dplyr)
library(ggplot2)

# Read data
loop_sizes <- read.table("loop_sizes_start_end.txt", 
                         header = TRUE, sep = "\t")
# Add loop ID for clarity
loop_sizes$loop_id <- seq_len(nrow(loop_sizes))
head(loop_sizes)

# Remove specific outliers
loop_sizes <- loop_sizes[loop_sizes$eupsc_size != 12900000, ]

# Pivot longer to get one row per species per loop
long_data <- pivot_longer(loop_sizes, cols = c(eupsc_size, sepof_size, octbi_size),
                          names_to = "species", values_to = "loop_size", values_drop_na = TRUE)

# Add genome size per species
genome_sizes <- c(eupsc_size = 5.1, sepof_size = 5.8, octbi_size = 2.7)
long_data$genome_size <- genome_sizes[long_data$species]

# Nicer species names for legend
long_data$species_clean <- recode(long_data$species,
                                  eupsc_size = "E. scolopes",
                                  sepof_size = "S. officinalis",
                                  octbi_size = "O. bimaculoides")

# Test, for each conserved loop, whether its size increases with genome size across the species where it's conserved.
# Group by loop_id and calculate correlation or slope per loop
per_loop_summary <- long_data %>%
  group_by(loop_id) %>%
  filter(n() >= 2) %>%  # Only loops conserved in at least 2 species
  summarise(
    rho = cor(genome_size, loop_size, method = "spearman"),
    slope = if (n() >= 2) coef(lm(loop_size ~ genome_size))[2] else NA_real_,
    n_species = n()
  )

# View the result
head(per_loop_summary)

# Summarise trends
mean(per_loop_summary$slope, na.rm = TRUE)    # Avg slope

# Tests if there is  a linear scaling effect — i.e., does loop size change by a consistent number of bp per 1 Gb genome size increase?
wilcox.test(per_loop_summary$slope, mu = 0, alternative = "two.sided")

# Conserved loop scatter plot with genome size

# Make sure loop_id is a factor for coloring
long_data$loop_id <- as.factor(long_data$loop_id)

# Italic species names + genome sizes
species_breaks <- c(2.7, 5.1, 5.8)
species_labels <- c(
  expression(italic("O. bimaculoides")~"(2.7~Gb)"),
  expression(italic("E. scolopes")~"(5.1~Gb)"),
  expression(italic("S. officinalis")~"(5.8~Gb)")
)

# Keep species order consistent
long_data$species_clean <- factor(
  long_data$species_clean,
  levels = c("O. bimaculoides", "E. scolopes", "S. officinalis")
)

# Distinct color palette
long_data$loop_id <- factor(long_data$loop_id)
ids <- levels(long_data$loop_id)
n   <- length(ids)

make_distinct_hcl <- function(n, c = 90, l = 60) {
  hues <- seq(0, 360, length.out = n + 1)[-1]
  grDevices::hcl(h = hues %% 360, c = c, l = l)
}

pal_distinct <- if (n <= 36) {
  make_distinct_hcl(n, c = 95, l = 60)
} else {
  lums <- rep(c(55, 70), length.out = n)
  hues <- seq(0, 360, length.out = n + 1)[-1]
  grDevices::hcl(h = hues %% 360, c = 100, l = lums)
}

# Randomize color order (so colours are easier to see)
set.seed(48) #Change if you want different colours
pal_random <- sample(pal_distinct, size = n, replace = FALSE)

# Plot
species_breaks <- c(2.7, 5.1, 5.8)
species_labels <- c(
  expression(italic("O. bimaculoides")~"(2.7~Gb)"),
  expression(italic("E. scolopes")~"(5.1~Gb)"),
  expression(italic("S. officinalis")~"(5.8~Gb)")
)

p <- ggplot(long_data, aes(x = genome_size, y = loop_size,
                           group = loop_id, color = loop_id)) +
  geom_line(size = 0.9, alpha = 0.9) +
  geom_point(size = 2.2, stroke = 0.3, alpha = 0.95) +
  scale_color_manual(values = pal_random, drop = FALSE) +
  scale_x_continuous(
    name = "Species (genome size)",
    breaks = species_breaks, labels = species_labels
  ) +
  labs(
    y = "Conserved loop size (bp)",
    title = "Conserved loop size vs. genome size",
    color = "Loop ID"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p

# Save to file
ggsave(
  filename = "conserved_loop_size_scatter.tiff",
  plot = p,
  device = "tiff",
  width = 8,
  height = 5,
  units = "in",
  dpi = 200
)

