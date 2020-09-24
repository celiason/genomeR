#' Save as single-line fasta
#' 
writefasta <- function(x, outpath, ext=".fa") {
	if (!is.null(names(unlist(x)))) {
		for (i in seq_along(x)) {
			cat(paste0(names(x[[i]]), "\n", x[[i]]), sep="\n", file=paste0(outpath, "/", names(x)[i], ext))
		}
	}
}

#' Abbreviate scientific names
#' 
#' @param x character of scientific names
#' 
#' @value output genome species abbreviates (e.g., "homSap" for Homo sapiens)
#' 
sppAbbrev <- function(x) {
	tmp <- str_split(x, " |_")
	paste0(tolower(substring(sapply(tmp, "[[", 1), 1, 3)),
		   toupper(substring(sapply(tmp, "[[", 2), 1, 1)),
		   substring(sapply(tmp, "[[", 2), 2, 3))
}
