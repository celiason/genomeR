# R function to plot gene structures output from exonerate

# source("https://bioconductor.org/biocLite.R")
# biocLite("Gviz")
library(Gviz)
library(Biostrings)

# plot genome ranges
seqs <- readDNAStringSet(paste0("data/genomes/todWin.fasta"))
names(seqs) <- str_extract(names(seqs), "^\\d+(?=\\s)")
genome(seqs) <- "todWin"
options(ucscChromosomeNames=FALSE)
sTrack <- SequenceTrack(seqs)
seqnames(sTrack)
chromosome(sTrack)
plotTracks(sTrack, chromosome="198", from=4429718, to=4429742)

# read files
files <- list.files("output/exonerate_models", "txt", full=T)
genes <- str_extract(files, "(?<=models\\/).*?(?=\\.txt)")
x <- lapply(files, readLines)

# extract GFF/gene structure data
gffs <- lapply(seq_along(x), function(i) {
	part1 <- x[[i]]
	start <- grep("START OF GFF DUMP", part1)
	end <- grep("END OF GFF DUMP", part1)
	part2 <- lapply(seq_along(start), function(z) { part1[start[z]:end[z]]})
	colnms <- strsplit(part2[[1]][grep("^\\d+", part2[[1]])[[1]]-2], " ")[[1]][-1]
	makedat <- function(x) {
		tf <- tempfile()
		writeLines(x[grep("^\\d+", x)], tf)
		res <- read.delim(tf, fill = TRUE, sep="\t", head=F)
		unlink(tf)
		names(res) <- colnms
		res
	}
	newres <- lapply(part2, makedat)
	newres
})

gffs[[5]]

names(gffs) <- genes

genes

library(seqinr)


# plot
res <- plotgenes(gffs, subset=3, cutoff="best")
res
# extract regions
extractSeq(gff=res, genome=seqs)


# loop - find/plot best gene models
pdf("figs/exonerate.pdf", width=8, height=3)
for (i in seq_along(genes)) {
	res <- plotgenes(gffs, subset=i, cutoff="best")
	myseqs <- extractSeq(gff=res, genome=seqs)
	write.fasta(myseqs, names=names(myseqs), file=paste0("output/exonerate_models/fasta/todWin_", genes[i], ".fa"))
}
dev.off()





write.fasta(myseqs, names=paste0(gr@seqnames, " GN=", gr$gene), file.out=paste0(spp, "_regions.fa"))

