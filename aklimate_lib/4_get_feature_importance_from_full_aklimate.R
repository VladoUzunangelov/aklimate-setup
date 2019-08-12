#!/usr/bin/env Rscript

## copyright (c) 2019 Vlado Uzunangelov

NUMBER_OF_CPUS_TO_USE = 28
GENERATE_FEATURE_IMPORTANCE_FILES = TRUE

message("Set some variables, load some libraries, source some paths.")

ncpus <- NUMBER_OF_CPUS_TO_USE
stopifnot(is.numeric(ncpus), ncpus > 0, ncpus%%1 == 0)
message("using ", ncpus, " CPUs")

modelsDir <- "./models"
reposDir <- "./repos"

target.dir <- getwd()

tasks <- 1:25
message("tasks: ", tasks)

library(doParallel)
library(foreach)
## library(doRNG)
registerDoParallel(cores = 15)

library(abind)
library(ranger)
library(dplyr)
library(caret)
library(ROCR)
library(pracma)
source(paste0(reposDir, "/tcga_scripts/utils.R"), chdir = TRUE)
source(paste0(reposDir, "/junkle/junkle-utils.R"), chdir = TRUE)
source(paste0(reposDir, "/junkle/junkle.R"), chdir = TRUE)
source(paste0(reposDir, "/Spicer/Spicer.R"), chdir = TRUE)
source(paste0(reposDir, "/Spicer/Spicer-classify.R"), chdir = TRUE)

library(FDb.InfiniumMethylation.hg19)

#######################################################

if (GENERATE_FEATURE_IMPORTANCE_FILES) {
  message("getting feature importance from full models")
  z1 <- foreach(i = iter(tasks)) %dopar% {
    load(paste0(modelsDir, "/", i, "_junkle_final_model.RData"))
    if (is.null(jklm)) {
      message("jklm object is null for ", i)
    }
    imps <- rank.features.jklm(jklm)
    if (is.null(imps)) {
      messge("imps is null for ", i)
    }
    write.df(data.frame(importance = imps), "features", paste0(modelsDir, "/",
      i, "_aklimate_multiclass_feature_importance.tab"))
  }
} else {
  message("skipping generating importance files")
}

#######################################################

message("DONE!")
