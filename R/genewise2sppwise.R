#' Convert gene-wise fasta files to species-wise fastas
#' 
#' Function converts gene-wise fasta files to species-wise fastas
#' (e.g., multiple sp. per one gene file)
#' Relies on faSomeRecords
#' Download at http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faSomeRecords
#' Put in /bin folder
#' Add execute permissions chmod a+x faSomeRecords.
#' 
#' @param files list of FASTA files
#' @param cores number of cores
#' @param contig name of contig to pull out
#' @param outpath path to output FASTAs
#' @param namesep unused?
#' 
genewise2sppwise <- function(files, cores=1, contig=NULL, outpath=".", namesep="|") {
	require(parallel)
	require(stringr)
	require(pbapply)
	if (!dir.exists("temp")) {
		dir.create("temp")
	}
	# Extract a specific contig:
	if (!is.null(contig)) {
		cat(contig, file="temp/ids.txt", sep="\n")
		# Extract contig fastas
		mclapply(files, mc.cores=cores, function(f) {
			# f <- files[1]
			nm <- gsub("\\.(fa|fasta)", paste0(".", contig, ".fa"), basename(f))
			system(paste0("faSomeRecords ", f, " temp/ids.txt temp/", nm))
		})
		subfiles <- list.files("temp", pattern=".fa", full=TRUE)
		nms <- substr(basename(subfiles), 1, 6)
		mclapply(seq_along(subfiles), mc.cores=cores, function(i) {
		  x <- readLines(subfiles[i])
		  x[[1]] <- gsub(paste0(">", contig), paste0(">", contig, "_", nms[i]), x[[1]])
		  x[[1]] <- gsub("\\s.*?$", "", x[[1]])
		  cat(x, file=subfiles[i], sep="\n")
		})
		# concatenate into single MSA file
		system(paste0("cat temp/*.fa > contig", contig, ".msa"))
		# cleanup
		system("rm temp/*")
	}
	seqs <- lapply(files, readLines)
	seqstarts <- lapply(seqs, grep, pattern=">")
	if (length(unique(sapply(seqstarts, length)))!=1) {
		stop("Not same number of genes")
	}
	nseq <- length(seqstarts[[1]])
	nsamp <- length(seqs)
	nms <- na.omit(str_extract(seqs[[1]], paste0("(?<=\\", namesep, ")(.*?)$")))
	# Now parse sequences by gene and output as new files
	for (i in 1:nseq) {
		for (j in 1:nsamp) {
			start <- seqstarts[[j]][i]
			if (i == length(seqstarts[[j]])) {
				end <- length(seqs[[j]])
			} else {
				end <- seqstarts[[j]][i+1] - 1
			}
			bits <- seqs[[j]][start : end]
			cat(bits, file = paste0(outpath, "/", nms[i], ".fa"), append=TRUE, sep="\n")
		}
	}
}
