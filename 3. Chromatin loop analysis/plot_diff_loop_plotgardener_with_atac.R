#Plotgardener code to plot differential loop for Figure 4 at stage 29, with ATAC-seq tracks.

rm(list = ls())

#plotgardner with juicebox files
#https://phanstiellab.github.io/plotgardener/articles/guides/reading_data_for_plotgardener.html

hicFile <- file.path("212493_intrachrom_rm_scaffolds.allValidPairs.hic")

# Load necessary packages
library(BSgenome)
library(BSgenomeForge)
library(devtools)
library("plotgardener")
library("plotgardenerData")
library("GenomicFeatures")
library("OrganismDbi")
library("AnnotationHub")
library(RColorBrewer)

#Make TxDb object from a gff/gtf for assembly.See https://rdrr.io/bioc/GenomicFeatures/man/makeTxDbFromGFF.html
eupsc_gff <- file.path("gmap.belcaid.cds.gff")

eupsc_txdb <- makeTxDbFromGFF(eupsc_gff,
                              format=c("auto"),
                              dataSource=NA,
                              organism= "Euprymna scolopes",
                              taxonomyId=NA,
                              circ_seqs=NULL,
                              chrominfo=NULL,
                              miRBaseBuild=NA,
                              metadata=NULL)


#Load orgdb that Oleg already made for GO analyses.
library(org.Eesc.eg.db)

# Create the genome assembly object
# See https://phanstiellab.github.io/plotgardener/reference/assembly.html

my_assembly <- assembly(
  Genome = "Lachesis_assembly.fasta",
  TxDb = eupsc_txdb,
  gene.id.column = "GENEID",
  display.column = "GENEID",
  OrgDb = org.Eesc.eg.db,
  BSgenome = "BSgenome.Eupsc.thea.eupscv2") # Provide the BSgenome object.


# Show the assembly and transcript db
print(my_assembly)
print(eupsc_txdb)

# Open a TIFF device to save the plot
tiff(
  filename = "diff_dev_loop_eupsc_50k_st29.tiff", 
  width = 4, 
  height = 3.5, 
  units = "in",    # Specify the units
  res = 300        # Resolution in DPI
)

## Create a page. Won't work without this grid, but you can remove it when you are finished.
pageCreate(width = 4, height = 2.5, default.units = "inches")

## Plot and place triangle Hi-C plot
hicPlot <- plotHicTriangle(
  data = hicFile, resolution = 50000,
  zrange = c(0, 5),
  chrom = "Lachesis_group4__56_contigs__length_174030884",
  chromstart = 60000000, chromend = 65000000,
  assembly = my_assembly,
  x = 2, y = 0.5, width = 3, height = 1.5,
  just = "top", default.units = "inches", palette = colorRampPalette(brewer.pal(n = 9, "Oranges")))

# Add "60Mb" on the left
plotText(
  label = "60Mb",
  x = 0.5, y = 2.03, just = c("left", "top"),
  fontsize = 8, fontface = "bold", col = "black",
  default.units = "inches"
)

# Add "Chromosome 5" in the middle
plotText(
  label = "Chromosome 5",
  x = 2, y = 2.03, just = c("center", "top"),  # Centered at x = 2.5
  fontsize = 8, fontface = "bold", col = "black",
  default.units = "inches"
)

# Add "65Mb" on the right
plotText(
  label = "65Mb",
  x = 3.5, y = 2.03, just = c("right", "top"),
  fontsize = 8, fontface = "bold", col = "black",
  default.units = "inches"
)

## Annotate heatmap legend
annoHeatmapLegend(
  plot = hicPlot, x = 3.5, y = 0.5,
  width = 0.13, height = 1.2,
  just = c("right", "top")
)

# Hide page guides
pageGuideHide()

# Prioritized gene order and highlights (include only preferred genes)
gene_pref <- c("cluster_1064.path1", "cluster_9293.path1", "cluster_3151.path1","cluster_11478.path1","cluster_4212.path1","cluster_18817.path1") #This only works with one gene when I input it as gene order..
gene_highlights <- data.frame(
  gene = c("cluster_1064.path1", "cluster_9293.path1", "cluster_3151.path1","cluster_11478.path1","cluster_4212.path1","cluster_18817.path1"),
  color = c("darkorange2", "darkorange2", "darkorange2","darkorange2", "darkorange2", "darkorange2"),
  stringsAsFactors = FALSE
)

## Plot gene track in same genomic region. Do fontsize = 0 to take off the labels.
genesPlot <- plotGenes(
  assembly = my_assembly, chrom = "Lachesis_group4__56_contigs__length_174030884", chromstart = 60000000, chromend = 65000000,
  x = 0.5, y = 2.1, width = 3, height = 0.5,
  just = c("left", "top"), default.units = "inches", geneHighlights = gene_highlights, fontsize = 0
)


#Add ATAC-seq data. Must be in BigWig format. 
atac_data <- import("97309_st29_chr5.bw")
atac_data_other <- import("97308_st25_chr5.bw")
atac_data_other2 <- import("97307_st20_chr5.bw")

# Pick a shared max (99.5th percentile across all three tracks in the region so they are scaled to the same height
qfun <- function(x, p=0.995) if (length(mcols(x)$score)) quantile(mcols(x)$score, p, na.rm=TRUE) else 0
ymax <- max(qfun(atac_data), qfun(atac_data_other), qfun(atac_data_other2))
yr <- c(0, ymax)

atacPlot3 <- plotSignal(
  data = atac_data_other2, chrom = "Lachesis_group4__56_contigs__length_174030884", chromstart = 60000000, chromend = 65000000,
  x = 0.5, y = 2.48, width = 3, height = 0.5, 
  just = c("left", "top"), default.units = "inches", range = yr,
  linecolor = "goldenrod", negData = FALSE)

atacPlot2 <- plotSignal(
  data = atac_data_other, chrom = "Lachesis_group4__56_contigs__length_174030884", chromstart = 60000000, chromend = 65000000,
  x = 0.5, y = 2.48, width = 3, height = 0.5, 
  just = c("left", "top"), default.units = "inches", range = yr,
  linecolor = "#BB5566", negData = FALSE)

atacPlot <- plotSignal(
  data = atac_data, chrom = "Lachesis_group4__56_contigs__length_174030884", chromstart = 60000000, chromend = 65000000,
  x = 0.5, y = 2.48, width = 3, height = 0.5, 
  just = c("left", "top"), default.units = "inches", range = yr,
  linecolor = "#004488", negData = FALSE)

# Close the TIFF device
dev.off()
