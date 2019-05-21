#!/bin/bash

#chrisw 20190520
#install various R packages for running AKLIMATE

R -e 'install.packages(c("foreach","doParallel","ranger","plyr","abind","ROCR","caret","proxy","purrr","pracma","fastmatch","devtools","mlr","e1071","igraph","circlize","RColorBrewer"))'
R -e 'source("https://bioconductor.org/biocLite.R")'
R -e 'BiocInstaller::biocLite(c("ComplexHeatmap"))'
R -e 'devtools::install_github("VladoUzunangelov/similarity")'
