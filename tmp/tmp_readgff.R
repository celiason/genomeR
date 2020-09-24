# Read GFFs in R??

if (!requireNamespace("BiocManager", quietly = TRUE) {
	install.packages("BiocManager")	
})

# BiocManager::install("GenomicFeatures")
# BiocManager::install("BiocParallel")

library(GenomicFeatures)

# GenomicFeatures::transcriptsByOverlaps

txdb <- makeTxDbFromGFF(file = "/Users/chadeliason/Desktop/todChl_rnd3.all.maker.protein2genome.gff", dataSource="testing", organism="Todiramphus chloris")

txdb

extractUpstreamSeqs(genome, txdb)

cdsBy(txdb, by = "gene")

gff <- read.delim(file = "/Users/chadeliason/Desktop/todChl_rnd3.all.maker.protein2genome.gff", skip=1)

dim(gff)

