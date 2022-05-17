# genomeR

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
