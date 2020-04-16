#!/bin/bash

#chrisw 20190520
#install various R packages for running AKLIMATE

# SPICER needs:
R -e 'install.packages("Rcpp", type = "source")'
R -e 'install.packages(c("plyr", "RcppArmadillo"))'
R -e 'devtools::install_github("VladoUzunangelov/SPICER")'
