# Load required library
library(ggplot2)

# Create data frame with line breaks in labels for better readability
data <- data.frame(
  species = factor(
    c(
      "E. scolopes loops",
      "S. officinalis loops",
      "O. bimaculoides loops"
    ),
    levels = c(
      "E. scolopes loops",
      "S. officinalis loops",
      "O. bimaculoides loops"
    )
  ),
  percentages = c(49.1, 57.9, 71.3)
)

#dark_blue_shades <- c(
 # "navy",        # classic navy blue
  #"darkblue",    # another shade of dark blue
  #"midnightblue", # a very dark, almost black blue
  #"slateblue",   # bluish grey with a dark tone
  #"royalblue4",  # darker variant of royal blue
  #"darkslateblue", # a mix of dark slate and blue
  #"indigo",      # a dark purplish-blue
  #"cornflowerblue", # a darker pastel blue, slightly muted
  #"steelblue4"   # darker variant of steelblue
# )


# Create the barplot
loop_chrom_stat_bar_pmax <- ggplot(data, aes(x = species, y = percentages, fill = species)) +
  geom_bar(stat = "identity", width = 0.6, alpha = 0.8) +  # Add alpha transparency
  scale_fill_manual(values = c("navy", "steelblue4", "cornflowerblue")) +
  ylim(0, 75) +
  labs(
    y = expression("Percentage of loops on different chromosome in "~italic("P. maxmimus")),
    x = NULL,
    title = " "
  ) +
  geom_text(aes(label = paste0(round(percentages, 0), "%")), vjust = -0.5, size = 5) +  # Adjusted text size
  theme_minimal() +
  theme(
    plot.margin = margin(10, 10, 10, 50),  # Adjust left margin
    plot.title = element_text(size = 18),  # Increase plot title size
    axis.text.x = element_text(face = "italic", size = 12, angle = 30, hjust = 1),  # Italic, smaller font, angled
    axis.text.y = element_text(size = 12),  # Adjust y-axis text size
    axis.title.y = element_text(size = 14),  # Adjust y-axis label size
    legend.position = "none"  # Remove legend
  )

# Print plot
print(loop_chrom_stat_bar_pmax)

ggsave("bar_loops_on_diff_chroms_pmax.tiff", loop_chrom_stat_bar_pmax, dpi = 200, width = 5, height = 7.5, units = "in")

