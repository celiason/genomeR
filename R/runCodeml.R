#' Function to run codeml within R
#' 
#' @param phy path to phylogenetic tree
#' @param model which model to use (see XX for details)
#' @param fasta path to FASTA codon alignment file
#' @param force whether to force analysis and mask early stop codons
#' @param outpath output path for results
#' @param silent whether to output run progress in console
#' 
#' @author Chad M. Eliason
#' 
#' testing
#' fasta=files[1]
#' model="M0"
#' phy="/home/FM/celiason/uce-alcedinidae/paml/kingtree_nolabels.phy"
#' 
runCodeml <- function(phy, model=c("M0", "branch-site", "free"), fasta, force=TRUE, outpath=NULL, silent=FALSE) {
	# fasta=f
	oldwd <- getwd()
	prefix <- gsub("\\..*?$", "", basename(fasta))
	if (!dir.exists("M0_results")) {
		dir.create("M0_results")
	}
	model <- match.arg(model)
	if (model=="branch-site") {
		x <- readLines("codeml-M2-BS-H1.ctl")
	}
	if (model=="free") {
		x <- readLines("codeml_free.ctl")
	}
	if (model=="M0") {
		x <- readLines("codeml-M0-BS.ctl")
		# ape::read.tree(phy)
		x <- stringr::str_replace(x, "treefile = .*", paste0("treefile = kingtree_nolabels.phy"))
		# warning("Ignoring `phy` argument. Using kingtree_nolabels.phy in working directory.")
		outpath <- paste0("M0_results/", prefix)
	}
	x <- stringr::str_replace(x, "treefile = .*", paste0("treefile = ", phy))
	x <- stringr::str_replace(x, "seqfile = .*", paste0("seqfile = ", fasta))
	outpath <- paste0("M0_results/", prefix)
	if (!dir.exists(outpath)) {
		dir.create(outpath)
	}
	setwd(outpath)
	cat(x, file="temp.ctl", sep="\n")
	# silent=FALSE
	system("codeml temp.ctl", ignore.stdout=silent)
	# system("yes", timeout=1) # automatically replace stop codons	
	res <- readLines("mlc")
	omega <- as.numeric(stringr::str_extract(grep("omega", res, value=T), "[\\d\\.]+"))
	setwd(oldwd)
	omega
}
