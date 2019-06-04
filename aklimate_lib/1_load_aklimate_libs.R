#!/usr/bin/env Rscript


## To be run on architeuthis in /home/ubuntu/projects/THYM_AKLIMATE

message("loading libraries")

ncpus <- 14
library(doParallel)
library(foreach)
## library(doRNG)
registerDoParallel(cores = ncpus)

library(abind)
library(ranger)
library(dplyr)
library(caret)
library(ROCR)
## library(mlr)
source("/aklimate_lib/tcga_scripts/utils.R", chdir = TRUE)
source("/aklimate_lib/junkle/junkle-utils.R", chdir = TRUE)
source("/aklimate_lib/junkle/junkle.R", chdir = TRUE)
source('/aklimate_lib/Spicer/Spicer.R',chdir=TRUE)
source('/aklimate_lib/Spicer/Spicer-classify.R',chdir=TRUE)

## need to install mySQL first sudo apt-get install libmariadbclient-dev then in R
## source('https://bioconductor.org/biocLite.R')
## biocLite('FDb.InfiniumMethylation.hg19')

library(FDb.InfiniumMethylation.hg19)
