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
#' @param namesep unused? e.g., "|" for ">todChl|Mc1r" format
#' 
genewise2sppwise <- function(files, cores=1, contig=NULL, outpath=".", group=FALSE, namesep=NULL, force=FALSE) {
	require(parallel)
	require(stringr)
	require(pbapply)
	# Extract a specific contig:
	if (!is.null(contig)) {
		if (!dir.exists("temp")) {
			dir.create("temp")
		}
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

	if (is.null(contig)) {
		seqs <- pblapply(files, readLines)
		seqstarts <- lapply(seqs, grep, pattern=">")
		if (!group & length(unique(sapply(seqstarts, length)))!=1) {
			stop("Not same number of genes. Consider grouping by name with `group`.")
		} #else if group {
		# 	str_extract(seqs[[1]], "(GMOD.*?$)")
		# 	grep("GMOD.*?", seqs[[1]], value=T)
		# }
		nseq <- length(seqstarts[[1]])
		nsamp <- length(seqs)
		# namesep="|"
		# seqs[[1]][1]
		nms <- gsub(">", "", seqs[[1]][seqstarts[[1]]])
		nms <- gsub("(-|\\.)", "_", nms)
		# nms <- na.omit(stringr::str_extract(seqs[[1]], paste0("(?<=\\", namesep, ")(.*?)$")))
		# Now parse sequences by gene and output as new files
		for (i in 1:nseq) {
			outfile <- paste0(outpath, "/", nms[i], ".fa")
			if (file.exists(outfile) & !force) {
				stop("Files exists, consider using force=TRUE?")
			} else if (force) {
				file.remove(outfile)
			}
			for (j in 1:nsamp) {
				start <- seqstarts[[j]][i]
				if (i == length(seqstarts[[j]])) {
					end <- length(seqs[[j]])
				} else {
					end <- seqstarts[[j]][i+1] - 1
				}
				spp <- substr(basename(files[j]), 1, 6)
				bits <- seqs[[j]][start : end]
				bits <- gsub(">.*?$", paste0(">", spp), bits)
				cat(bits, file=outfile, append=TRUE, sep="\n")
			}
		}
	}
}
