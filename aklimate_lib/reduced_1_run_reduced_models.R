#!/usr/bin/env Rscript

## copyright (c) 2019 Vlado Uzunangelov

message("Set some variables, load some libraries, source some paths.")

ncpus <- NUMBER_OF_CPUS_TO_USE
stopifnot(is.numeric(ncpus), ncpus > 0, ncpus%%1 == 0)
message("using ", ncpus, " CPUs")

dataDir <- "/home/ubuntu/remote_dirs/courtyard/pstore_migration/home_dir/mkl/data/tmp_awg/LIHC"
message("dataDir: ", dataDir)

target.dir <- getwd()
message("target.dir: ", target.dir)

tasks <- 1:25
message("tasks: ", tasks)

library(doParallel)
library(foreach)
## library(doRNG)
registerDoParallel(cores = ncpus)

library(abind)
library(ranger)
library(dplyr)
library(caret)
library(ROCR)
library(pracma)
source("repos/tcga_scripts/utils.R", chdir = TRUE)
source("repos/junkle/junkle-utils.R", chdir = TRUE)
source("repos/junkle/junkle.R", chdir = TRUE)
source("repos/Spicer/Spicer.R", chdir = TRUE)
source("repos/Spicer/Spicer-classify.R", chdir = TRUE)

library(FDb.InfiniumMethylation.hg19)


############################################################################

# This section copied from '2_prepare_data.R'.

if (exists("dat") && exists("pathways") && exists("labels")) {
  message("data already exists. skip data loading.")
} else {
  message("data objects do not exist. need to read and preprocess.")

  QUANTIZE_NUMERIC_DATA = FALSE

  message("load sample data")

  ## MIR is also there - about 800 MIRs right now probably best to not include them
  ## - I have a scheme to indluce them by putting a MIR in every pathway that has
  ## one of its targets, but it didn't work well when I tested it
  nomir = c("MUTA", "CNVR", "METH", "GEXP")
  suffs <- list(nomir)

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


  message("filter out low expression features")
  # Theo has named all gene expression features with 'N:GEXP:'

  sampleIDs <- rownames(dat)
  exp_matrix <- dat[sampleIDs, grepl(paste0("N:GEXP:"), colnames(dat))]

  message(paste("number of features in input expression matrix", length(colnames(exp_matrix))))

  exp_means <- colMeans(exp_matrix)
  exp_sds <- apply(exp_matrix, 2, sd)

  ## exclude genes that are in the bottom quartile by mean or sd
  exp_features_upper75_by_mean = names(exp_means)[exp_means > quantile(exp_means,
    probs = 0.25)]
  exp_features_upper75_by_sd = names(exp_sds)[exp_sds > quantile(exp_sds, probs = 0.25)]
  exp_all_features <- colnames(exp_matrix)
  exp_keep_features <- intersect(exp_features_upper75_by_mean, exp_features_upper75_by_sd)
  exp_drop_features <- setdiff(exp_all_features, exp_keep_features)

  full_keep_features <- setdiff(colnames(dat), exp_drop_features)

  message(paste("number of expression features to drop", length(exp_drop_features)))

  exp_matrix <- NULL

  message(paste("number of features in full data matrix", length(colnames(dat))))
  dat <- dat[sampleIDs, full_keep_features]
  message(paste("number of features in new data matrix", length(colnames(dat))))


  if (QUANTIZE_NUMERIC_DATA) {
    message("quantize numeric datatypes")

    dat <- foreach(datatype_str = iter(nomir[!grepl("MUTA|CNVR", nomir)]), .combine = cbind) %do%
      {
        message(paste("quantize ", datatype_str))
        quantize.data(dat[sampleIDs, grepl(paste0("N:", datatype_str, ":"),
          colnames(dat))], nbr = 5, idx = sampleIDs)
      }
  }

  # aklimate expects a dataframe, not a matrix
  dat <- as.data.frame(dat)

  message("load pathways")

  homeDir <- "./p_store_files"
  workDir <- "./models/"

  # p1 <- readSetList(paste0(homeDir, '/pathcomm_pathways_cleaned')) p1 <-
  # p1[!grepl('[[:digit:]]_SMPDB$|Z_SMPDB$', names(p1))]
  p1 <- readSetList(paste0(homeDir, "/pathcomm_pathways_cleaned_non_redundant_names.tsv"))
  p2 <- readSetList(paste0(homeDir, "/genomic_position_sets.listt"))
  p3 <- readSetList(paste0(homeDir, "/genesigdb_human.tab"))
  p4 <- readSetList(paste0(homeDir, "/msigdb_c2_c5_no_c2_cp.tab"))

  pathways <- c(p1, p2, p3, p4)
  max_size_of_pathways <- 1000
  pathways <- pathways[sapply(pathways, length) < max_size_of_pathways]


  message("sanitize pathway names")

  ## probably a good idea to sanitize the name a bit
  names(pathways) <- gsub("[^\\w\\s]", "_", names(pathways), perl = TRUE)

  message(paste("number of pathways to use", length(names(pathways))))

  message("load labels")

  labels <- as.matrix(read.delim("./labels.tsv", check.names = F, stringsAsFactors = F,
    header = F, row.names = 1))
  labels <- factor(labels[, 1])
  lbls <- labels


  message("load CV fold splits")

  splits <- as.matrix(read.delim("./data/cv_folds.tsv", check.names = F, stringsAsFactors = F,
    header = TRUE, row.names = 1))
  splits <- splits[, -1]


  message("set tasks and seeds")

  tasks <- colnames(splits)
  tasks <- tasks[1:25]

  # you can change the seeds - this is just to match initial runs
  seeds <- 11 * (1:length(tasks))
  names(seeds) <- tasks

  message("FINISHED - set tasks and seeds")

}
#############################

# labels <- as.matrix(read.delim(paste0(dataDir, '/iclust_membership'),
# check.names = F, stringsAsFactors = F, header = T, row.names = 1))
# load(paste0(dataDir,
# '/20180731_combined_features_binned_5_filtered_mean_var.RData'))
# load(paste0(dataDir,
# '/path_comm_pathways_methylation_mirna_extended_plus_positional_sets.RData'))
# dat <- dat.comb.binned labels <- labels[rownames(dat), ] pathways <-
# path.comm.pathways.positional.sets pathways <-
# pathways[!grepl('[[:digit:]]_SMPDB$|Z_SMPDB$', names(pathways))] groups <-
# c('exp', 'methylation', 'mutation', 'cnv', 'cn_regions', 'lncrna') new_folds <-
# readSetList(paste0(dataDir, '/cv_folds.listt')) new_universe <-
# scan(paste0(dataDir, '/all_samples'), what = character()) lbls <-
# factor(labels) names(lbls) <- names(labels) lbls <- lbls[new_universe] dat <-
# dat[new_universe, ]

#######################################################


#
z1 <- foreach(i = iter(tasks)) %do% {
  # load(paste0(target.dir, '/split_', i, '_aklimate_multiclass_model.RData'))
  task <- colnames(splits)[i]
  message("z1 starting ", i)
  model_file_path <- paste0(target.dir, "/models/", task, "_junkle_final_model.RData")

  message("model_file_path: ", model_file_path)
  load(model_file_path)

  message("z1 load completed ", i)
  imps <- rank.features.jklm(jklm)

  message("z1 rank.features.jklm completed ", i)
  feature_importance_file_path <- paste0(target.dir, "/models/", task, "_aklimate_multiclass_feature_importance.tab")

  message("feature_importance_file_path: ", feature_importance_file_path)
  write.df(data.frame(importance = imps), "features", feature_importance_file_path)

  message("z1 output written ", i)
}

cutoffs <- c(1, 5, 10, 50, 100, 1000)
reps.list <- lapply(1:25, function(x) seq(5 * (x - 1) + 1, 5 * x, 1))

stopifnot(FALSE)

acc.reduced <- foreach(i = iter(reps.list), .combine = rbind) %dopar% {
  foreach(k = iter(cutoffs), .combine = c) %do% {
    preds <- foreach(j = iter(i), .combine = c) %do% {
      idx.train <- new_folds[[j]]
      idx.test <- setdiff(new_universe, idx.train)



      imps <- read.delim(paste0(target.dir, "/split_", j, "_aklimate_multiclass_feature_importance.tab"),
        header = TRUE, row.names = 1)
      dat <- mlr::createDummyFeatures(dat)
      rf <- ranger(data = cbind(data.frame(labels = lbls[idx.train]), dat[idx.train,
        rownames(imps)[1:k], drop = FALSE]), dependent.variable.name = "labels",
        always.split.variables = NULL, classification = TRUE, sample.fraction = 0.5,
        num.trees = 3000, mtry = ceiling(k/5), min.node.size = 1, case.weights = NULL,
        num.threads = 3, probability = TRUE, respect.unordered.factors = FALSE,
        importance = "none", write.forest = TRUE, keep.inbag = TRUE, replace = FALSE)
      save(rf, file = paste0(target.dir, "/split_", j, "_cutoff_", k, "_rf_reduced_model.RData"))
      rf.preds <- predict(rf, dat[idx.test, rownames(imps)[1:k]])$predictions
      rownames(rf.preds) <- idx.test
      save(rf.preds, file = paste0(target.dir, "/split_", j, "_cutoff_", k,
        "_rf_reduced_model_predictions.RData"))

      confM <- caret::confusionMatrix(factor(apply(rf.preds[idx.test, ], 1,
        which.max), levels = levels(lbls)), lbls[idx.test])
      save(confM, file = paste0(target.dir, "/split_", j, "_cutoff_", k, "_rf_reduced_model_stats.RData"))

      mean(unname(confM$byClass[, "Balanced Accuracy"]))
      ## mean(unname(confM$overall['Accuracy']))
    }

    mean(preds)

  }
}

colnames(acc.reduced) <- cutoffs
write.df(data.frame(acc.reduced, check.names = FALSE), "repetition", paste0(target.dir,
  "/balanced_accuracy_reduced_feature_sets.tab"))

######################################################
acc.reduced <- foreach(i = iter(reps.list), .combine = rbind) %dopar% {
  foreach(k = iter(cutoffs), .combine = c) %do% {
    preds <- foreach(j = iter(i), .combine = c) %do% {
      load(paste0(target.dir, "/split_", j, "_cutoff_", k, "_rf_reduced_model_stats.RData"))
      mean(unname(confM$overall["Accuracy"]))
    }

    mean(preds)

  }
}

colnames(acc.reduced) <- cutoffs
write.df(data.frame(acc.reduced, check.names = FALSE), "repetition", paste0(target.dir,
  "/accuracy_reduced_feature_sets.tab"))

####################################################################### same as above, but with feature selection proportions derived from aklimate
####################################################################### weights
cutoffs <- c(5, 10, 20, 50, 100, 200, 500, 1000, 1500)
reps.list <- lapply(1:25, function(x) seq(5 * (x - 1) + 1, 5 * x, 1))


acc.reduced.akl <- foreach(i = iter(reps.list), .combine = rbind) %dopar% {
  foreach(k = iter(cutoffs), .combine = c) %do% {
    preds <- foreach(j = iter(i), .combine = c) %do% {
      idx.train <- new_folds[[j]]
      idx.test <- setdiff(new_universe, idx.train)



      imps <- read.delim(paste0(target.dir, "/split_", j, "_aklimate_multiclass_feature_importance.tab"),
        header = TRUE, row.names = 1)
      dat <- mlr::createDummyFeatures(dat)

      rf <- ranger(data = cbind(data.frame(labels = lbls[idx.train]), dat[idx.train,
        rownames(imps)[1:k], drop = FALSE]), dependent.variable.name = "labels",
        always.split.variables = NULL, classification = TRUE, sample.fraction = 0.5,
        num.trees = 3000, mtry = ceiling(k/5), min.node.size = 1, case.weights = NULL,
        split.select.weights = imps[1:k, 1]/sum(imps[1:k, 1]), num.threads = 3,
        probability = TRUE, respect.unordered.factors = FALSE, importance = "none",
        write.forest = TRUE, keep.inbag = TRUE, replace = FALSE)
      save(rf, file = paste0(target.dir, "/split_", j, "_cutoff_", k, "_rf_reduced_model_aklimate_weighting.RData"))
      rf.preds <- predict(rf, dat[idx.test, rownames(imps)[1:k]])$predictions
      rownames(rf.preds) <- idx.test
      save(rf.preds, file = paste0(target.dir, "/split_", j, "_cutoff_", k,
        "_rf_reduced_model_aklimate_weighting_predictions.RData"))

      confM <- caret::confusionMatrix(factor(apply(rf.preds[idx.test, ], 1,
        which.max), levels = levels(lbls)), lbls[idx.test])
      save(confM, file = paste0(target.dir, "/split_", j, "_cutoff_", k, "_rf_reduced_model_aklimate_weighting_stats.RData"))

      mean(unname(confM$byClass[, "Balanced Accuracy"]))
      ## mean(unname(confM$overall['Accuracy']))
    }

    mean(preds)

  }
}

colnames(acc.reduced.akl) <- cutoffs
write.df(data.frame(acc.reduced.akl, check.names = FALSE), "repetition", paste0(target.dir,
  "/balanced_accuracy_reduced_feature_sets_aklimate_weighting.tab"))





###########################################
overall.acc <- foreach(i = iter(reps.list), .combine = c) %dopar% {
  scores <- foreach(k = iter(i), .combine = c) %do% {

    load(paste0(target.dir, "/split_", k, "_predictions.RData"))
    idx.train <- new_folds[[k]]
    idx.test <- setdiff(new_universe, idx.train)

    confM <- caret::confusionMatrix(factor(jklm.preds, levels = levels(lbls)),
      lbls[idx.test])

    mean(unname(confM$byClass[, "Balanced Accuracy"]))
    ## mean(unname(confM$overall['Accuracy']))

  }
  mean(scores)
}

write.df(data.frame(overall = overall.acc), "repetition", paste0(target.dir, "/balanced_accuracy_full_aklimate.tab"))
## write.df(data.frame(overall=overall.acc),'repetition',paste0(target.dir,'/accuracy_full_aklimate.tab'))

confM.list <- foreach(i = iter(reps.list)) %dopar% {
  scores <- foreach(k = iter(i)) %do% {

    load(paste0(target.dir, "/split_", k, "_predictions.RData"))
    idx.train <- new_folds[[k]]
    idx.test <- setdiff(new_universe, idx.train)

    confM <- caret::confusionMatrix(factor(jklm.preds, levels = levels(lbls)),
      lbls[idx.test])

    confM$table
  }

  confM.avg <- Reduce("+", scores)
  confM.avg <- as.data.frame.matrix(confM.avg)
  rownames(confM.avg) <- colnames(confM.avg) <- paste0("iClust_", colnames(confM.avg))
  confM.avg

}

save(confM.list, file = paste0(target.dir, "/confusion_matrix_list.RData"))

confM.avg <- Reduce("+", confM.list)/length(confM.list)
confM.avg <- as.data.frame.matrix(confM.avg)

write.df(confM.avg, "", paste0(target.dir, "/confusion_matrix_average.tab"))


######################################################################## 3

props.list <- foreach(i = iter(reps.list)) %dopar% {
  props <- foreach(j = iter(i)) %do% {
    imps <- read.delim(paste0(target.dir, "/split_", j, "_aklimate_multiclass_feature_importance.tab"),
      header = TRUE, row.names = 1)
    imps <- setNames(imps[, 1], rownames(imps))
    coffs <- c(cutoffs, length(imps))
    types.prop <- foreach(k = iter(coffs), .combine = cbind) %do% {
      out <- rank.importance.type(groups, imps[1:k], intron = "_")
      out[groups]
    }
    colnames(types.prop) <- c(as.character(coffs[-length(coffs)]), "full")
    types.prop
  }
  Reduce("+", props)/length(props)
}
avg.props <- Reduce("+", props.list)/length(props.list)
avg.props <- avg.props[order(avg.props[, "full"], decreasing = TRUE), ]

write.df(data.frame(avg.props, check.names = FALSE), "data_type", paste0(target.dir,
  "/data_type_contributions.tab"))

######################################################### 3

t1 <- foreach(i = iter(reps.list), .combine = cbind) %dopar% {
  t2 <- foreach(j = iter(i), .combine = c) %do% {
    idx.train <- new_folds[[j]]
    idx.test <- setdiff(new_universe, idx.train)
    load(paste0(target.dir, "/split_", j, "_predictions.RData"))
    t3 <- as.numeric(jklm.preds)
    names(t3) <- names(jklm.preds)
    t3
  }
  t2 <- t2[new_universe]
  t2
}

write.df(t1, "sample", paste0(target.dir, "/aklimate_cv_individual_calls_334_samples.tab"))

t1 <- foreach(i = iter(reps.list), .combine = cbind) %dopar% {
  t2 <- foreach(j = iter(i), .combine = c) %do% {
    idx.train <- new_folds[[j]]
    idx.test <- setdiff(new_universe, idx.train)
    load(paste0(target.dir, "/split_", j, "_predictions.RData"))
    jklm.preds.bin <- as.numeric(as.numeric(jklm.preds[idx.test]) == as.numeric(lbls[idx.test]))
    names(jklm.preds.bin) <- idx.test
    jklm.preds.bin
  }
  t2 <- t2[new_universe]
  t2
}

stats <- cbind(ncol(t1), rowMeans(t1))
colnames(stats) <- c("count", "proportion")
write.df(stats, "sample", paste0(target.dir, "/aklimate_cv_per_sample_accuracy_334_samples.tab"))

######################################################################################

cf <- 20
## cf <- 50
t1 <- foreach(i = iter(reps.list), .combine = cbind) %dopar% {
  t2 <- foreach(j = iter(i), .combine = c) %do% {
    load(paste0(target.dir, "/split_", j, "_cutoff_", cf, "_rf_reduced_model_predictions.RData"))
    t3 <- apply(rf.preds, 1, function(x) as.numeric(colnames(rf.preds)[which.max(x)]))
    names(t3) <- rownames(rf.preds)
    t3
  }
  t2 <- t2[new_universe]
  t2
}

write.df(data.frame(t1), "sample", paste0(target.dir, "/rf_reduced_model_cutoff_",
  cf, "_individual_calls_334_samples.tab"))

t1 <- foreach(i = iter(reps.list), .combine = cbind) %dopar% {
  t2 <- foreach(j = iter(i), .combine = c) %do% {
    idx.train <- new_folds[[j]]
    idx.test <- setdiff(new_universe, idx.train)
    load(paste0(target.dir, "/split_", j, "_cutoff_", cf, "_rf_reduced_model_predictions.RData"))
    t3 <- apply(rf.preds, 1, function(x) as.numeric(colnames(rf.preds)[which.max(x)]))
    names(t3) <- rownames(rf.preds)
    t3 <- as.numeric(t3[idx.test] == as.numeric(lbls[idx.test]))
    names(t3) <- idx.test
    t3
  }
  t2 <- t2[new_universe]
  t2
}
stats <- cbind(ncol(t1), rowMeans(t1))
colnames(stats) <- c("count", "proportion")
write.df(stats, "sample", paste0(target.dir, "/rf_reduced_model_cutoff_", cf, "_per_sample_accuracy_334_samples.tab"))

####################################################################
