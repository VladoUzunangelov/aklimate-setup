# A toy example for running AKLIMATE on a binary classification task

# vuzunangelov & chrisw 20210210

message("==> begin load libraries")
NUMBER_OF_CPUS_TO_USE <- 1
ncpus <- NUMBER_OF_CPUS_TO_USE
message("using ", ncpus, " CPUs")


message("==> use AKLIMATE package")

library(aklimate)

message("==> setup parallelization")
# Internally, AKLIMATE uses %dopar% for parallelization
library(doParallel)
stopifnot(ncpus > 0, ncpus <= detectCores())
registerDoParallel(cores = ncpus)
number_parWorkers <- getDoParWorkers()
message(paste0("==> number of parWorkers: ", number_parWorkers))

message("==> caret for confusion matrix")
library(caret)

message("==> utility functions")
readSetList <- function(file, delim = "\t") {
  l <- readLines(file)
  l.fields <- strsplit(l, delim)
  r <- lapply(l.fields, function(x) as.vector(x[-1]))
  names(r) <- sapply(l.fields, "[[", 1)
  return(r)
}

write.df <- function(df, row.names.id = "", out.file) {
  output <- cbind(rownames(df), df)
  colnames(output)[1] <- row.names.id
  write.table(output, file = out.file, quote = FALSE, append = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)
}

message("==> finish load libraries")
message("==> begin load data")


message("==> setup dirs")
workDir <- getwd()

task.dir <- paste(workDir, "example_results", sep = "/")
if (!dir.exists(task.dir)) dir.create(task.dir)

num_cv_folds <- 4
tasks <- 1:num_cv_folds
# tasks <- c(1)


message("==> pathways")
pathways_file <- paste0(workDir, "/pathways.tsv")
pathways <- readSetList(pathways_file)

message("==> labels")
labels_file <- paste0(workDir, "/labels.tsv")
labels <- as.matrix(read.delim(labels_file, check.names = F, stringsAsFactors = F,
  header = T, row.names = 1))

message("==> features")
sample_data_file <- paste0(workDir, "/sample_data.tsv")
sample_data <- as.data.frame(read.delim(sample_data_file, check.names = F, stringsAsFactors = F,
  header = T, row.names = 1))

message("==> splits")
splits <- createFolds(labels[, 1], k = num_cv_folds, list = TRUE, returnTrain = TRUE)

message("==> finish load data")
message("==> begin running AKLIMATE")


suffix <- "_75_25_aklimate_model.RData"


groups <- list(exp = c("exp"))

dat <- sample_data


message("==> tasks")
message(tasks)

###################################################################
worker.f <- function(tasks) {
  res <- foreach(i = iter(tasks)) %do% {
    set.seed(11 * i, kind = "L'Ecuyer-CMRG")

    message("==> setting parameters. Some settings used here sacrifice classification performance for speed.")
    message("description of parameters at https://github.com/VladoUzunangelov/aklimate")

    idx.train <- rownames(labels)[splits[[i]]]
    idx.test <- setdiff(rownames(labels), idx.train)

    training_data <- dat[idx.train, ]
    datatypes <- groups
    training_lables <- labels[idx.train, 1]
    feature_sets <- pathways
    global_features_names <- NULL

    rf_params <- list(ttype = "binary", bin.perf = "bacc", regression.q = 0.05,
      importance = "impurity", min.nfeat = 15, ntree = 2000, sample.frac = 0.5,
      replace = FALSE, weights = NULL, oob.cv = data.frame(min.node.prop = 0.01,
        mtry.prop = 0.25, ntree = 500))

    aklimate_params <- list(topn = 5, subsetCV = TRUE, lamb = c(-15, 0), cvlen = 20,
      celnet = c(0.01, 0.002), type = "probability")

    store_kernels <- FALSE
    verbose <- TRUE

    message("==> train model")

    ta <- Sys.time()

    aklimate_model <- aklimate(training_data, datatypes, training_lables, feature_sets,
      global_features_names, rf_params, aklimate_params, store_kernels, verbose)

    tb <- Sys.time()

    t_train <- difftime(tb, ta, units = "secs")
    message(paste0("training time: ", t_train))

    save(aklimate_model, file = paste0(task.dir, "/split_", i, suffix))
    message("==> saved model")

    message("==> test model")
    tc <- Sys.time()
    aklimate_model.preds <- predict(aklimate_model, dat, pathways, NULL, FALSE)$preds

    td <- Sys.time()

    t_test <- difftime(td, tc, units = "secs")
    message(paste0("testing time: ", t_test))

    message(paste0("train + test time: ", t_train + t_test))

    factor.test <- as.factor(labels[idx.test, 1])
    factor.preds <- aklimate_model.preds[, 2]
    factor.preds[aklimate_model.preds[, 2] < 0.5] <- colnames(aklimate_model.preds)[1]
    factor.preds[aklimate_model.preds[, 2] > 0.5] <- colnames(aklimate_model.preds)[2]
    factor.preds <- factor(factor.preds, levels = levels(factor.test))
    confM <- confusionMatrix(factor.preds, factor.test)
    rocr.pred <- ROCR::prediction(aklimate_model.preds[, 2], labels[idx.test,
      1])
    auc <- ROCR::performance(rocr.pred, "auc")@y.values[[1]]

    save(aklimate_model.preds, file = paste0(task.dir, "/split_", i, "_predictions.RData"))
    save(confM, rocr.pred, auc, file = paste0(task.dir, "/split_", i, "_stats.RData"))

    message("==> saved test results")

    gc()

    return(list(auc = auc))

  }


  return(res)
}

message("==> train/test models")
t0 <- Sys.time()
results <- lapply(tasks, worker.f)
t1 <- Sys.time()
message("==> got results object")
dt <- difftime(t1, t0, units = "secs")
message(dt)
dt

results <- unlist(results, recursive = FALSE)

##########################
pred_stats <- sapply(results, function(x) x[[1]])
names(pred_stats) <- as.character(tasks)

write.df(data.frame(auc = pred_stats), "split", paste0(task.dir, "/splits_auc.tab"))




######################################################

pred_stats <- foreach(i = tasks, .combine = c) %do% {
  load(paste0(task.dir, "/split_", i, "_stats.RData"))
  auc
}
names(pred_stats) <- tasks

message("==> finish running AKLIMATE")
