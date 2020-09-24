# Function to extract consensus sequences from gff file
bam2fasta <- function(gff, bamid=NULL, ref="King_genome/kingfisher.fasta") {
	for (j in unique(gff$gene)) {
		gffsub <- subset(gff, gene==j)
		raw <- list()
		for (i in 1:nrow(gffsub)) {
			region <- paste0(gffsub$seqname[i], ":", gffsub$start[i], "-", gffsub$end[i])
			if (!is.null(bamid)) {
				cat("Processing gene", i)
				system(paste0("bcftools mpileup -q 30 -Q 30 -d 100 -r ", region, " -Ou -f King_genome/kingfisher.fasta ", bamid, ".bam | bcftools call -mv -Oz -o calls.vcf.gz"))
				# normalize indels
				system(paste0("bcftools norm -f King_genome/kingfisher.fasta -m +any -Oz calls.vcf.gz -o calls.norm.vcf.gz"))
				# filter adjacent indels within 5bp
				system(paste0("bcftools index -f calls.norm.vcf.gz"))
				raw[[i]] <- system(paste0("samtools faidx King_genome/kingfisher.fasta ", region, " | bcftools consensus -H 1 calls.norm.vcf.gz"), intern=TRUE)
				# raw[[i]] <- system(paste0("samtools faidx King_genome/kingfisher.fasta ", region, " | bcftools consensus -i TYPE='snp' -H 1 calls.norm.vcf.gz"), intern=TRUE)
			} else {
				raw[[i]] <- system(paste0("samtools faidx King_genome/kingfisher.fasta ", region), intern=TRUE)
			}
			raw[[i]][1] <- paste0(">", bamid, "_", gffsub$gene[i], "_", gsub(">","",raw[[i]][1]))
		}
		rawc <- paste0(do.call(paste0, lapply(raw, "[", -1)), collapse="")
		# Format as fasta (60 wide columns)
		rawc <- substring(rawc, seq(1, nchar(rawc), 60), seq(60, nchar(rawc), 60))
		# Combine single fasta files to multifasta
		if (!is.null(bamid)) {
			cat(rawc, file=paste0(bamid,"_taste_genes.fa"), append=TRUE, sep="\n")
			# Add sequence/gene name
			rawc <- c(paste0(">", bamid, "_", j), rawc)
			system("rm calls*")
		} else {
			rawc <- c(paste0(">reference_", j), rawc)
			cat(rawc, file="reference_taste_genes.fa", append=TRUE, sep="\n")
		}
	}
}