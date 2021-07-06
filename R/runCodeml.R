#' Function to run codeml within R
#' 
#' @param phy path to phylogenetic tree
#' @param model which model to use (see XX for details)
#' @param fasta path to FASTA codon alignment file
#' @param force whether to force analysis and mask early stop codons
#' @param outpath output path for results
#' @param silent whether to output run progress in console
#' @param M0path path for M0 results to use for branch lengths
#' 
#' @author Chad M. Eliason
#' 
#' testing
#' fasta=files[1]
#' model="M0"
#' phy="/home/FM/celiason/uce-alcedinidae/paml/kingtree_nolabels.phy"
#' 
runCodeml <- function(fasta, model=c("M0", "M1a", "M2a", "M7", "M8", "BS", "free", "MA", "MAnull"), phy, force=FALSE, fixedbl=FALSE, outpath=NULL, silent=FALSE, M0path=NULL) {
	# fasta=f
	require(stringr)
	require(treeio)
	oldwd <- getwd()
	prefix <- gsub("\\..[^\\.]*$", "", basename(fasta))
	
	model <- match.arg(model)
	# if (model=="BS") {
	# 	ctl <- readLines("codeml-M2-BS-H1.ctl")
	# }
	# if (model=="free") {
	# 	ctl <- readLines("free.ctl")
	# }

	# TODO: create MA, MAnull ctl files
	ctl <- readLines(paste0(model, ".ctl"))

	# For model M0 we are using phylogeny without any nodes labeled
	if (model == "M0") {
		ctl <- str_replace(ctl, "treefile = .*", paste0("treefile = kingtree_nolabels.phy"))
		warning("Ignoring `phy` argument. Using kingtree_nolabels.phy in working directory.")
	}

	if (is.null(M0path)) {
		M0path <- paste0("M0_results/", prefix)	
	} else {
		M0path <- paste0(M0path, "/", prefix)
	}
	
	# Setup outut path for files
	if (is.null(outpath)) {
		outpath <- paste0(model, "_results/", prefix)		
	} else {
		outpath <- paste0(outpath, "/", prefix)
	}

	# For models M1a, M2a, and free-ratio we are using branch lengths estimated under M0 model:
	if (model %in% c("M1a", "M2a", "M7", "M8", "free", "MA", "MAnull")) {
		raw <- readLines(paste0(M0path, "/mlc"))
		# raw <- readLines(M0_mlc)
		textphy <- raw[grep("tree length", raw)[1]+4]
	}

	# if (dir.exists(paste0(model, "_results/", prefix))) {
	# 	stop(paste0("Path ", paste0(model, "_results/", prefix), " already exists, consider deleting?"))
	# } else {
	# 	dir.create(paste0(model, "_results/", prefix))
	# }

	if (!dir.exists(outpath)) {
		dir.create(outpath, recursive=TRUE)
	} else if (dir.exists(outpath) & !force) {
		stop(paste0("Path ", outpath, " already exists, consider deleting?"))
	}
	
	setwd(outpath)

	if (fixedbl) {
		cat(textphy, file="fixed.phy")
		if (model %in% c("MA", "MAnull")) {
			# modify the M0 trees and add annotations for nodes (e.g., "#1")
			tr <- treeio::read.tree("fixed.phy")
			tr$node.label <- NULL
			# plot(tr, show.node.label=TRUE)
			# ape::nodelabels()
			tr <- treeio::label_branch_paml(tr, 44, "#1")
			tr <- treeio::label_branch_paml(tr, 56, "#1")
			tr <- treeio::label_branch_paml(tr, 59, "#1")
			tr <- treeio::label_branch_paml(tr, 51, "#1")
			# plot(tr, show.node.label=TRUE)
			treeio::write.tree(tr, file="fixed.phy")
		}
		ctl <- str_replace(ctl, "treefile = .*$", "treefile = fixed.phy")
	} else {
		ctl <- str_replace(ctl, "treefile = .*", paste0("treefile = ", phy))	
	}

	ctl <- str_replace(ctl, "seqfile = .*", paste0("seqfile = ", fasta))
	cat(ctl, file="temp.ctl", sep="\n")

	# Run
	system("codeml temp.ctl", ignore.stdout=silent)

	setwd(oldwd)

}

parseCodeml <- function(respath) {
	list.files(respath, pattern="mlc")
	res <- readLines("mlc")
	omega <- as.numeric(stringr::str_extract(grep("omega", res, value=T), "[\\d\\.]+"))
	setwd(oldwd)
	omega
}
