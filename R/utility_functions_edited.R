# Function to subset the seqinr list
# Extracts sequences from a genome based on coordinates
extractSeq <- function(gff, genome) {
# gff <- res
	# if (all(gff$strand=="-")) {
	# 	gff <- arrange(gff, end)
	# }
	# if (all(gff$strand=="+")) {
	# 	gff <- arrange(gff, start)
	# }
	nm <- as.character(gff$seqname)  # what scaffold it's on
	id <- match(nm, names(genome))
	res <- lapply(1:nrow(gff), function(i) {
		x <- gff[i, ]
		ss <- id[i]
		if (as.character(x$strand) == '-') {
      		reverseComplement(genome[[ss]][x$start:x$end])
    	} else {
      	genome[[ss]][x$start:x$end]
    	}    	
	})
	names(res) <- paste0(gff$seqname, "_", gff$start, "_", gff$end)
	res
}

getgene <- function(gene, uids, outpath=NULL) {
	x <- paste0('txid', uids, '[Organism:noexp] AND ', gene, '[Gene] AND ("1"[SLEN] : "8000"[SLEN])')
	x <- lapply(x, entrez_search, db='protein', retmax=1)
	x <- x[!as.numeric(sapply(x, "[[", "count"))==0]
	ids <- sapply(x, "[[", "ids")
	x <- sapply(1:length(ids), function(x) {entrez_fetch(db="nuccore", id = ids[[x]], rettype="fasta")})
	x <- gsub('\n\n', '\n', x)
	if (!is.null(outpath)) {
		outpath <- ifelse(grepl("\\/$", outpath), outpath, paste0(outpath, "/"))
		cat(x, file=paste0(outpath, "/", gene, ".fa"), sep="\n")
	} else {
		cat(x, sep="\n")
	}
	
}


# x = list of gff data.frames
getscores <- function(x) {
	as.numeric(as.character(sapply(lapply(x, "[[", "score"), "[", 1)))
}


# file=
# cutoff= score cutoff to use
plotgenes <- function(gff, subset, cutoff="best") {
# gff=gffs
# subset=18
	require(Gviz)
	genename <- names(gff[subset])
	gff <- gff[[subset]]
	scores <- getscores(gff)
	if (cutoff=="best") {
		keep <- which.max(scores)
	} else if (is.numeric(cutoff)) {
		keep <- scores > cutoff	
		if (all(scores < cutoff)) {
			stop("Adjust cutoff, too low")
		}
	}
	gff <- gff[keep]
	chrom <- rep(seq_along(gff), sapply(gff, nrow))
	gff <- do.call(rbind, gff)
	# gff$chromosome <- chrom
	# length(gffs[[16]])
	# gff$feature <- plyr::revalue(gff$feature, c('intron'='utr'))
	gff$exon <- NA
	gff$transcript <- gff$seqname
	gff$symbol <- genename
	gff$gene <- chrom
	gff$lengths <- abs(gff$end-gff$start)
	gff <- subset(gff, feature == "exon")
	grt <- GeneRegionTrack(gff, name=genename)
	plotTracks(grt, showId=TRUE, collapseTranscripts=FALSE)
	invisible(gff)
}