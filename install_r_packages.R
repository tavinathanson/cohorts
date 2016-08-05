# commands to install needed R packages

source("http://bioconductor.org/biocLite.R")
biocLite("BSgenome")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
biocLite("GenomeInfoDb")
biocLite("VariantAnnotation")

cran_mirror <- "http://cran.us.r-project.org"

install.packages('plyr', repos=cran_mirror, dependencies=TRUE)
install.packages('Rcpp', repos=cran_mirror, dependencies=TRUE)
install.packages('testthat', repos=cran_mirror, dependencies=TRUE)
install.packages('knitr', repos=cran_mirror, dependencies=TRUE)
install.packages('reshape2', repos=cran_mirror, dependencies=TRUE)
install.packages('deconstructSigs', repos=cran_mirror, dependencies=TRUE)