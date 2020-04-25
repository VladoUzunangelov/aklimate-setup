#!/bin/bash

#chrisw 200200424
#install various R packages for running AKLIMATE

# SPICER needs:
R -e 'install.packages("Rcpp", type = "source")'
R -e 'install.packages(c("plyr", "RcppArmadillo"))'

# aklimate should pull in spicer
# R -e 'devtools::install_github("VladoUzunangelov/SPICER")'

R -e 'devtools::install_github("VladoUzunangelov/aklimate")'
