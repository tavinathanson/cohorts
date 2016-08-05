# commands to install needed R packages

bio_source <- "https://bioconductor.org/biocLite.R"
source(bio_source)
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