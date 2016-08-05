# commands to install needed R packages

source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
biocLite("GenomeInfoDb")
biocLite("VariantAnnotation")

install.packages('plyr')
install.packages('Rcpp')
install.packages('testthat')
install.packages('knitr')
install.packages('reshape2', dependencies=TRUE)

cran_mirror <- "http://cran.us.r-project.org"
install.packages('deconstructSigs', repos=cran_mirror, dependencies=TRUE)