# Genome analysis in R

This is the development version of the R package genomeR

I wrote these functions to streamline comparative genomics projects I've been working on. Use at your own risk, as this package is currently under development.

You can install the package using the devtools R package:

```r
install.packages("devtools")
devtools::install_github("celiason/genomeR")
library(genomeR)
```

Then, for example, to align reads in file "myreads.fastq.gz" to reference "myref.fa" run the following:

```r
alignReads(ref="myref.fa", reads="myreads.fastq", cores=48, ram=150, suffix="run1")
```

Note: the directory structure should be such that there is a folder called "genomes" with 6-letter species abbreviations as subfolders (e.g., "anaPla" for the mallard, _Anas platyrhynchos_). The code above will run a QC analysis using the `fastp` program (must be installed) and output a sorted BAM alignment file with the following naming convention: READS-to-REF.bam (e.g., "anaPla-to-galGal.bam" for mallard aligned to the chicken genome _Gallus gallus_).
