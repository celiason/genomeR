#' Function to align raw reads to reference genome
#' 
#' @param ref path to reference genome (genome will be indexed if not already)
#' @param reads path to raw reads file (PE)
#' @param cores how many cores to use for BWA MEM
#' @param ram how much RAM (GB) to allot to BAM file output
#' 
alignReads <- function(ref, reads, cores=48, ram=150) {
	# Paths to programs
	bwa="/home/FM/celiason/bwa/bwa"
	fastp="/home/FM/celiason/fastp"
	# Setup
	oldwd <- getwd()
	# id <- 
	id <- gsub(".fa", "", basename(ref))
	from <- substr(basename(reads), 1, 6)
	to <- substr(basename(ref), 1, 6)
	prefix <- paste0(from, "-to-", to)
	# Index reference genome
	setwd(paste0("genomes/", from))
	if (!file.exists(paste0(id, ".fa.bwt"))) {
		system(paste0(bwa, " index ", basename(ref)))
	}
	setwd(oldwd)
	# Align raw reads to reference genome
	if (file.exists(paste0("alignments/", prefix, ".bam"))) {
		stop("BAM file already exists!")
	}
	system(paste0(fastp, " -i ", reads, " --interleaved_in --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --stdout -h alignments/", prefix, ".html | sed -E 's/^((@|\\+)SRR[^.]+\\.[^.]+)\\.(1|2)/\\1/' | ", bwa, " mem -p -t ", cores, " ", ref, " - | samtools view -Sb - | samtools sort -m ", ram, "G > alignments/", prefix, ".bam"))
}