# Load required library
library(ggplot2)

# Create data frame with line breaks in labels for better readability
data <- data.frame(
  species = factor(
    c(
      "E. scolopes loops, different\nO. bimaculoides chromosomes",
      "S. officinalis loops, different\nO. bimaculoides chromosomes",
      "O. bimaculoides loops, different\nE. scolopes chromosomes",
      "O. bimaculoides loops, different\nS. officinalis chromosomes"
    ),
    levels = c(
      "E. scolopes loops, different\nO. bimaculoides chromosomes",
      "S. officinalis loops, different\nO. bimaculoides chromosomes",
      "O. bimaculoides loops, different\nE. scolopes chromosomes",
      "O. bimaculoides loops, different\nS. officinalis chromosomes"
    )
  ),
  percentages = c(22.22, 24.93, 52.42, 56.90)
)

#Other colour options "plum3", "#785EF0", "#A38DF5", "#A95AA1"

# Create the barplot
loop_chrom_stat_bar <- ggplot(data, aes(x = species, y = percentages, fill = species)) +
  geom_bar(stat = "identity", width = 0.6, alpha = 0.8) +  # Add alpha transparency
  scale_fill_manual(values = c("#117733", "#44AA99", "#882255", "plum3")) +
  ylim(0, 75) +
  labs(
    y = "Percentage of loops",
    x = NULL,
    title = " "
  ) +
  geom_text(aes(label = paste0(round(percentages, 0), "%")), vjust = -0.5, size = 5) +  # Adjusted text size
  theme_minimal() +
  theme(
    plot.margin = margin(10, 10, 10, 50),  # Adjust left margin
    plot.title = element_text(size = 18),  # Increase plot title size
    axis.text.x = element_text(face = "italic", size = 10, angle = 30, hjust = 1),  # Italic, smaller font, angled
    axis.text.y = element_text(size = 12),  # Adjust y-axis text size
    axis.title.y = element_text(size = 14),  # Adjust y-axis label size
    legend.position = "none"  # Remove legend
  )

# Print plot
print( loop_chrom_stat_bar)
ggsave("bar_loops_on_diff_chroms.tiff", loop_chrom_stat_bar, dpi = 200, width = 7.5, height = 7.5, units = "in")
