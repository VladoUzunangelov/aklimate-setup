#!/usr/bin/env Rscript
# vuzunangelov & chrisw 20190606

## Run the commands in this script after those in 1_load_aklimate_libs.R.
## example directory tree:
# .
# ├── 1_load_aklimate_libs.R
# ├── 2_run_aklimate.R
# ├── data
# │   ├── combined_matrix.tsv
# │   └── cv_folds.tsv
# ├── datatypes.tsv
# ├── files.txt
# ├── labels.tsv
# ├── models
# ├── p_store_files
# │   ├── genesigdb_human.tab
# │   ├── genomic_position_sets.listt
# │   ├── msigdb_c2_c5_no_c2_cp.tab
# │   └── pathcomm_pathways_cleaned
# ├── repos
# │   ├── junkle
# │   │   ├── junkle.R
# │   │   └── junkle-utils.R
# │   ├── Spicer
# │   │   ├── experimental
# │   │   │   ├── Isomap.R
# │   │   │   ├── SPARKLE.R
# │   │   │   ├── Spicer-red.R
# │   │   │   └── ToDo
# │   │   ├── kernels
# │   │   │   ├── CumulativeRBF.R
# │   │   │   ├── CustomRBF.R
# │   │   │   └── Kernels.R
# │   │   ├── LICENSE
# │   │   ├── prank
# │   │   │   └── pRank.R
# │   │   ├── README.md
# │   │   ├── Spicer-classify.R
# │   │   ├── Spicer-funcs.cpp
# │   │   ├── Spicer-funcs.R
# │   │   ├── Spicer.R
# │   │   ├── ToDo
# │   │   └── utils.R
# │   └── tcga_scripts
# │       └── utils.R
# ├── run_aklimate.mak
# └── samples.tsv

message("load sample data")

## MIR is also there - about 800 MIRs right now probably best to not include them
## - I have a scheme to indluce them by putting a MIR in every pathway that has
## one of its targets, but it didn't work well when I tested it

suffs <- list(nomir = c("MUTA", "CNVR", "METH", "GEXP"))

dat <- read.delim("./data/combined_matrix.tsv", header = T, row.names = 1, check.names = FALSE)
dat <- dat[, -1]


message("meth mapper")

## load('../data/meth_mappers.RData') move this to a makefile or something - take
## a couple of minutes to compute and should only be done once
hm450 <- get450k()
meth.mapper <- getNearestGene(hm450)

updated.names <- foreach(i = iter(colnames(dat)), .combine = c) %dopar% {
  if (grepl("METH", i, ignore.case = TRUE)) {
    cg <- strsplit(i, ":")[[1]][4]
    gsub(paste0(cg, "::"), paste0(cg, ":", meth.mapper[cg, 4], ":"), i)
  } else {
    i
  }
}

colnames(dat) <- updated.names



message("load pathways")

homeDir <- "./p_store_files"
workDir <- "./models/"


p1 <- readSetList(paste0(homeDir, "/pathcomm_pathways_cleaned"))
p1 <- p1[!grepl("[[:digit:]]_SMPDB$|Z_SMPDB$", names(p1))]
p2 <- readSetList(paste0(homeDir, "/genomic_position_sets.listt"))
p3 <- readSetList(paste0(homeDir, "/genesigdb_human.tab"))
p4 <- readSetList(paste0(homeDir, "/msigdb_c2_c5_no_c2_cp.tab"))

pathways <- c(p1, p2, p3, p4)
pathways <- pathways[sapply(pathways, length) < 1000]


message("sanitize pathway names")

## probably a good idea to sanitize the name a bit
names(pathways) <- gsub("[^\\w\\s]", "_", names(pathways), perl = TRUE)



message("load labels")

labels <- as.matrix(read.delim("./labels.tsv", check.names = F, stringsAsFactors = F,
  header = F, row.names = 1))
labels <- factor(labels[, 1])


message("load CV fold splits")

splits <- as.matrix(read.delim("./data/cv_folds.tsv", check.names = F,
  stringsAsFactors = F, header = TRUE, row.names = 1))
splits <- splits[, -1]


message("set tasks and seeds")

tasks <- colnames(splits)

# you can change the seeds - this is just to match initial runs
seeds <- 11 * (1:length(tasks))
names(seeds) <- tasks



message("define worker.f")

worker.f <- function(tasks) {

  res <- foreach(i = iter(tasks)) %do% {

    set.seed(seeds[i], kind = "L'Ecuyer-CMRG")

    idx.train <- rownames(splits)[splits[, i] == 0]
    idx.test <- setdiff(rownames(splits), idx.train)

    jklm <- junkle(dat[idx.train, ], suffs, labels[idx.train], pathways, NULL,
      list(ttype = "multiclass", bin.perf = c("bacc"), importance = "permutation",
        min.nfeat = 15, ntree = 1000, sample.frac = 0.5, replace = FALSE,
        weights = NULL, oob.cv = data.frame(min.node.prop = 0.01, mtry.prop = 0.25,
          ntree = 500)), list(topn = 5, subsetCV = TRUE, lamb = c(-8, 3),
        cvlen = 200, type = "probability"), FALSE, TRUE)

    save(jklm, file = paste0(workDir, "/", i, "_junkle_final_model.RData"))

    jklm.preds <- predict.junkle(jklm, dat[c(idx.train, idx.test), ], pathways,
      NULL, FALSE)$preds
    jklm.preds <- apply(jklm.preds, 1, function(x) colnames(jklm.preds)[which.max(x)])

    confM <- caret::confusionMatrix(factor(jklm.preds, levels = levels(labels)),
      labels[idx.test])
    bacc <- mean(unname(confM$byClass[, "Balanced Accuracy"]))



    save(jklm.preds, confM, bacc, file = paste0(workDir, "/", i, "_junkle_final_model_stats_preds.RData"))

    return(bacc)


  }

  return(res)
}



message("call worker on tasks")

# EDIT THIS BEFORE RUNNING !
run_tasks <- tasks[1:25]
stopifnot(length(run_tasks) > 0)
message(paste(length(run_tasks), " tasks to run"))
message(run_tasks)


# stopifnot(FALSE)

results <- worker.f(run_tasks)
names(results) <- run_tasks

results <- unlist(results)

names(results) <- as.character(run_tasks)



message("write results")

write.df(data.frame(bacc = results), "split", paste0(workDir, "/splits_balanced_accuracy.tab"))
