#!/bin/bash

#chrisw 20190520
#install various R packages for running AKLIMATE

R -e 'install.packages(c("foreach","doParallel","ranger","plyr","abind","ROCR","caret","proxy","purrr","pracma","fastmatch","devtools","mlr","e1071","igraph","circlize","RColorBrewer","RMySQL"))'
R -e 'source("https://bioconductor.org/biocLite.R")'
R -e 'BiocInstaller::biocLite(c("ComplexHeatmap","FDb.InfiniumMethylation.hg19"))'
R -e 'devtools::install_github("VladoUzunangelov/similarity")'

#TODO do we also need "/home/ubuntu/repos/tcga_scripts/utils.R" ?
