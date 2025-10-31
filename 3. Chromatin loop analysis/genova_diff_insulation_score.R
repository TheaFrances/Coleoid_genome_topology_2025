library(GENOVA)
library(ggplot2)

rm(list = ls())

# Load contacts
eupsc_20_50kb <- load_contacts(
  signal_path = '/320995_50000_iced.matrix',
  indices_path = '320995_50000_abs.bed',
  sample_name = "Stage 20",
  balancing = "T",
  colour = "goldenrod"
)

eupsc_25_50kb <- load_contacts(
  signal_path = '212492_50000_iced.matrix',
  indices_path = '212492_50000_abs.bed',
  sample_name = "Stage 25",
  balancing = "T",
  colour = "#BB5566"
)

eupsc_29_50kb <- load_contacts(
  signal_path = '212493_50000_iced.matrix', 
  indices_path = '212493_50000_abs.bed',
  sample_name = "Stage 29",
  balancing = "T",
  colour = "#004488"
)


# Generate insulation score
options(datatable.allow.cartesian=TRUE)
eupsc_20v25v29_50kb_ins <- insulation_score(
  list(eupsc_20_50kb, eupsc_25_50kb, eupsc_29_50kb),
  window = 13 # 10 corresponds to a 500kb window at 50kb resolution
)


# Create the initial plot
plot_insulation <- visualise(eupsc_20v25v29_50kb_ins,
                             chr = 'Lachesis_group4__56_contigs__length_174030884', start = 60e6, end = 65e6,
                             contrast = 1)

# Add annotations with transparency to the highlighted regions
plot_with_annotations <- plot_insulation +
  geom_rect(data = data.frame(xmin = c(61644484, 61746148, 61772744, 63297532, 63342407, 63342440),
                              xmax = c(61724980, 61774884, 61826740, 63339246, 63435852, 63435670),
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "darkorange", alpha = 0.5, inherit.aes = FALSE) + # Apply transparency only to the highlighted regions
  labs(title = "Insulation score with highlighted genes",
       x = "Position (bp)",
       y = "Insulation score") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  xlim(c(start = 60e6, end = 65e6))

# Save the plot with annotations
ggsave("eupsc_20vs25v29_50kb_insulation_dif_650kbwindow_diffloopdev2.tiff",
       plot = plot_with_annotations, units = "in", dpi = 200, width = 8, height = 6)

