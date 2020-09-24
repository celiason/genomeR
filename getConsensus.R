#' Function to get consensus sequences from BAM alignments
#' 
#' @param ref path to reference genome
#' @param bam path to BAM alignment file
#' @param call whether to call variants
#' @param index whether to index the VCF variant file output
#' @param cons whether to generate consensus sequence
#' 
#' Testing zone:
#' list.files("genomes/actHom")
#' ref="genomes/actHom/actHom-to-todChl.consensus.fa"
#' bam="alignments/actHom-to-actHom.bam"
#' 
getConsensus <- function(ref, bam, call=TRUE, index=TRUE, cons=TRUE, suffix=NULL) {
	nm <- stringr::str_extract(basename(bam), ".*?(?=\\.bam)")
	outpath <- paste0("genomes/", substr(nm, 1, 6))
	if (!dir.exists(outpath)) {
		dir.create(outpath)
	}
	# call variants (filter out indels with -I - Taylor Hains suggested doing this to maintain reference coordinates)
	if (call) {
		if (file.exists(paste0(outpath, "/", nm, ".calls.vcf.gz"))) {
			stop("VCF file already exists")
		}
		system(paste0("bcftools mpileup -I -Ou -f ", ref, " ", bam, " | bcftools call -mv -Oz -o ", outpath, "/", nm, ".calls.vcf.gz"))
	}
	# index
	if (index) {
		system(paste0("bcftools index ", outpath, "/", nm, ".calls.vcf.gz"))
	}
	# normalize indels
	# system(paste0("bcftools norm -f ", ref, " ", nm, ".calls.vcf.gz -Ob -o ", nm, ".calls.norm.bcf"))
	# filter adjacent indels within 5bp
	# system(paste0("bcftools filter --IndelGap 5 calls.norm.bcf -Ob -o calls.norm.flt-indels.bcf")
	if (cons) {
		if (file.exists(paste0(outpath, "/", nm, ".consensus.fa"))) {
			stop("Consensus FASTA already exists, try a different suffix")
		}
		system(paste0("cat ", ref, " | bcftools consensus ", outpath, "/", nm, ".calls.vcf.gz > ", outpath, "/", nm, ".consensus.fa"))
	}
	# This assumes we have already made the calls, normalized indels and filtered. There is another page which goes deeper and is devoted just to this, but in brief, the variant calling command in its simplest form is:
}
