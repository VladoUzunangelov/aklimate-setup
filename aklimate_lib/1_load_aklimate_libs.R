#!/usr/bin/env Rscript
# vuzunangelov & chrisw 20190606

## check 2_run_aklimate.R for example directory structure

message("loading libraries")

ncpus <- NUMBER_OF_CPUS_TO_USE
stopifnot(is.numeric(ncpus), ncpus > 0, ncpus%%1 == 0)
message("using ", ncpus, " CPUs")

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
aklimateWorkDir <- "/data"
source(paste0(aklimateWorkDir, "/repos/tcga_scripts/utils.R"), chdir = TRUE)
source(paste0(aklimateWorkDir, "/repos/junkle/junkle-utils.R"), chdir = TRUE)
source(paste0(aklimateWorkDir, "/repos/junkle/junkle.R"), chdir = TRUE)
source(paste0(aklimateWorkDir, "/repos/Spicer/Spicer.R"), chdir = TRUE)
source(paste0(aklimateWorkDir, "/repos/Spicer/Spicer-classify.R"), chdir = TRUE)

## need to install mySQL first sudo apt-get install libmariadbclient-dev then in R
## source('https://bioconductor.org/biocLite.R')
## biocLite('FDb.InfiniumMethylation.hg19')

library(FDb.InfiniumMethylation.hg19)
