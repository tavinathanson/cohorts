# commands to install needed R packages

source("http://bioconductor.org/biocLite.R")
biocLite("BSgenome", dependencies=TRUE)
biocLite("BSgenome.Hsapiens.UCSC.hg19", dependencies=TRUE)
biocLite("GenomeInfoDb", dependencies=TRUE)
biocLite("VariantAnnotation", dependencies=TRUE)

cran_mirror <- "http://cran.cnr.berkeley.edu/"

install.packages('plyr', repos=cran_mirror, dependencies=TRUE)
install.packages('Rcpp', repos=cran_mirror, dependencies=TRUE)
install.packages('testthat', repos=cran_mirror, dependencies=TRUE)
install.packages('knitr', repos=cran_mirror, dependencies=TRUE)
install.packages('reshape2', repos=cran_mirror, dependencies=TRUE)
install.packages('deconstructSigs', repos=cran_mirror, dependencies=TRUE)