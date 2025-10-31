# GO analysis for stage 29 E. scolopes loop genes

rm(list = ls())

library(dplyr)
library(ggplot2)
library(scales)  # For percentage formatting


# Make TxDb object from a gff/gtf for assembly.See https://rdrr.io/bioc/GenomicFeatures/man/makeTxDbFromGFF.html
eupsc_gff <- file.path("/gmap.belcaid.cds.gff")

eupsc_txdb <- makeTxDbFromGFF(eupsc_gff,
                              format=c("auto"),
                              dataSource=NA,
                              organism= "Euprymna scolopes",
                              taxonomyId=NA,
                              circ_seqs=NULL,
                              chrominfo=NULL,
                              miRBaseBuild=NA,
                              metadata=NULL)


#Load orgdb
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

# Take genes of interest and test their enrichment against whole set of genes according to different GO classifications (molecular functions, biological process etc)
# Import the go terms as a 'library' into R
library(clusterProfiler)
library(org.Eesc.eg.db)
library(ggplot2)
library(GO.db)
go2name=select(GO.db, keys=keys(GO.db), columns=columns(GO.db)[1:2])

# Your gene list you want to test for enrichment:
file="tot_loops_eupsc_29_50k+100k.tsv.genes_rm_dups_list.txt"

genes=read.table(file)
gs=unique(as.character(genes[,1]))

# Background gene list with GO terms made with createTerm2Gene.pl
allgenes=read.table("term2gene") 
ag=unique(as.character(allgenes[,2])) 


for (ont in c("CC","MF","BP")) {
  print (ont)
  
  ego <- enrichGO(gene         = gs,
                  universe      = ag,
                  keyType="GENENAME",
                  OrgDb         = org.Eesc.eg.db,
                  ont           = ont,
                  pAdjustMethod = "BH",
                  
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  write.table(ego,paste(file,"_enrich.",ont,".go",sep=''))
  
  if (nrow(ego)>0) {
    goplot(ego)
    ggsave(paste(file, "_enrichUnstable.",ont,".pdf",sep=''))
  }
  
}

for (ont in c("CC","MF","BP")) {
  print(ont)
  
  ego <- enrichGO(gene         = gs,
                  universe      = ag,
                  keyType       = "GENENAME",
                  OrgDb         = org.Eesc.eg.db,
                  ont           = ont,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  
  write.table(ego, paste(file, "_enrich.", ont, ".go", sep=''))
  
  if (nrow(ego) > 0) {
    p <- dotplot(ego, showCategory = 10) +  # `showCategory` controls how many GO terms are displayed
      theme(
        text = element_text(size = 8),         # General text size
        axis.text = element_text(size = 6),    # Axis tick labels
        axis.title = element_text(size = 8),   # Axis titles
        legend.text = element_text(size = 6),  # Legend text
        legend.title = element_text(size = 8), # Legend title
        plot.title = element_text(size = 10)   # Plot title
      ) +
      guides(size = guide_legend(override.aes = list(size = 3)))  # Makes legend points smaller
    
    ggsave(paste(file, "_enrichUnstable.", ont, ".pdf", sep=''), plot = p, width = 6, height = 6, dpi = 200)
  }
}




