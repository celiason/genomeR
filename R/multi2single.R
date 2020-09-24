#' Convert multi-line to single-line fastas
#' 
#' @param fastapath path to FASTA files
#' @param outpath output path
#' 
#' @export
#' 
multi2single <- function(fastapath, outpath=NULL) {
	files <- list.files(fastapath, pattern=".fa$", full=TRUE)
	seqs <- lapply(files, readLines)
	for (i in seq_along(seqs)) {
		# i=2
		firstlines <- grep(">", seqs[[i]])
		lastlines <- c(grep(">", seqs[[i]])[-1] - 1, length(seqs[[i]]))
		newseq <- lapply(seq_along(firstlines), function(x) {
			c(seqs[[i]][firstlines[x]], paste0(seqs[[i]][(firstlines[x]+1) : lastlines[x]], collapse=""))
		})
		newseq <- unlist(newseq)
		nms <- newseq[seq(1, length(newseq)-1, by=2)]
		seqs[[i]] <- setNames(as.list(newseq[seq(2, length(newseq), by=2)]), nms)
	}
	setNames(seqs, basename(files))
}
