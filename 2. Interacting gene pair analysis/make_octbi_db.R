#Make O. bimaculoides annotation library for R

library(AnnotationForge)

go=read.table("gid.go",header=T,sep="\t")
gene=read.table("gname.go",header=T,sep="\t")

makeOrgPackage( gene_info=gene, go=go,
                version="0.1",
                maintainer="octbi <test@test.com>",
                author="octbi <test@test.com>",
                outputDir = ".",
                tax_id="0",
                genus="octbi",
                species="octbi",
                goTable="go")


install.packages("./org.Ooctbi.eg.db", repos = NULL, type = "source")
