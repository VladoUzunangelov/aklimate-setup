#!/bin/bash

#chrisw 20190520
#install various R packages for running AKLIMATE

# ROCR needs "gplots", which needs "caTools".
# caTools_1.17.1.2 can be installed for R_3.4.4 with:
# require(devtools); install_version("caTools", version = "1.17.1.2")

R -e 'install.packages("devtools"); require(devtools); install_version("caTools", version = "1.17.1.2")'

R -e 'install.packages(c("foreach","doParallel","ranger","plyr","abind","ROCR","caret","proxy","purrr","pracma","fastmatch","devtools","mlr","e1071","igraph","circlize","RColorBrewer","RMySQL"))'
R -e 'source("https://bioconductor.org/biocLite.R")'
R -e 'BiocInstaller::biocLite(c("ComplexHeatmap","FDb.InfiniumMethylation.hg19"))'
R -e 'devtools::install_github("VladoUzunangelov/similarity")'

#TODO do we also need "/home/ubuntu/repos/tcga_scripts/utils.R" ?
