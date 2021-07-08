#!/usr/bin/env Rscript
# vuzunangelov & chrisw 20200521

LOAD_DATA <- !(exists("dat_raw") && exists("labels") && exists("pathways") && exists("cv_feature_whitelists") &&
  exists("tasks"))

TRAIN_TEST <- !(exists("bacc_results"))

LOAD_DATA <- TRUE
TRAIN_TEST <- TRUE


ncpus <- NUMBER_OF_CPUS_TO_USE
stopifnot(is.numeric(ncpus), ncpus > 0, ncpus%%1 == 0)
message("using ", ncpus, " CPUs")

cohortDir <- getwd()
featureSetsDir <- paste0(cohortDir, "/p_store_files")
modelsDir <- paste0(cohortDir, "/models")
dataDir <- paste0(cohortDir, "/data")

homeDir <- featureSetsDir



############################################################################


# sample_data_file <- paste0(dataDir, '/combined_matrix.tsv')

get_sample_data <- function(filename) {
  message(paste0("getting sample data from ", filename))

  dat_raw <- read.delim(filename, header = T, row.names = 1, check.names = FALSE)
  feature_data <- dat_raw[, c(-1)]

  sample_labels <- dat_raw[, 1]
  names(sample_labels) <- rownames(dat_raw)

  out <- list(sample_labels, feature_data)
  # get sample_labels with out[[1]] get feature_data with out[[2]]
}


############################################################################


detect_classification_type <- function(sample_labels) {

  num_labels <- length(levels(sample_labels))
  if (num_labels == 2) {
    classification_type = "binary"
  } else if (num_labels > 2) {
    classification_type = "multiclass"
  } else {
    message(paste0("ERROR: detected ", num_labels, " labels"))
    stopifnot(FALSE)
  }

  return(classification_type)
}


############################################################################


restrict_sample_data_to_feature_whitelist <- function(full_sample_data, cv_features_lists,
  cv_name) {

  dat_cols <- colnames(full_sample_data)
  # cv_name <- 'R1:F1'
  feature_whitelist <- cv_features_lists[cv_name][[1]]
  features_int <- intersect(dat_cols, feature_whitelist)
  dat_white <- full_sample_data[, features_int]
  message(paste0("restricted sample data\tsamples:", dim(dat_white)[1], "\tfeatures:",
    dim(dat_white)[2]))

  # return(dat_white)
  return(features_int)
}


############################################################################


meth_map_sample_data <- function(sample_data, meth_mapper) {
  message("renaming METH features to associate with gene symbol")
  updated.names <- foreach(i = iter(colnames(sample_data)), .combine = c) %dopar%
    {
      if (grepl("METH", i, ignore.case = TRUE)) {
        cg <- strsplit(i, ":")[[1]][4]
        gsub(paste0(cg, "::"), paste0(cg, ":", meth_mapper[cg, 4], ":"),
          i)
      } else {
        i
      }
    }

  colnames(sample_data) <- updated.names
}


############################################################################


# pathways_file <- paste0(featureSetsDir, '/collected_pathways_name_mapped.tsv')

load_pathways <- function(filename) {

  message("load pathways")

  # using one file of collected pathways.  The names of these pathways have been
  # mapped to simple names.  This was done to avoid the continued problem of name
  # collisions and the like.  collected_pathways_name_mapped.tsv
  pathways <- readSetList(pathways_file)

  max_size_of_pathways <- 1000
  pathways <- pathways[sapply(pathways, length) < max_size_of_pathways]

  return(pathways)
}


############################################################################


# cv_folds_file <- paste0(dataDir, '/cv_folds.tsv')

load_cv_splits <- function(filename) {
  message(paste0("loading cv splits from ", filename))
  splits <- as.matrix(read.delim(filename, check.names = F, stringsAsFactors = F,
    header = TRUE, row.names = 1))
  splits <- splits[, -1]

  return(splits)
}


############################################################################

if (LOAD_DATA) {
  cv_folds_file <- paste0(dataDir, "/cv_folds.tsv")
  splits <- load_cv_splits(cv_folds_file)

  message(paste0("loaded splits\tsamples:", dim(splits)[1], "\tsplits:", dim(splits)[2]))


  pathways_file <- paste0(featureSetsDir, "/collected_pathways_name_mapped.tsv")
  pathways <- load_pathways(pathways_file)

  message(paste0("num pathways loaded:", length(pathways)))

  message("load sample data")

  ## MIR is also there - about 800 MIRs right now probably best to not include them
  ## - I have a scheme to indluce them by putting a MIR in every pathway that has
  ## one of its targets, but it didn't work well when I tested it
  nomir = c("MUTA", "CNVR", "METH", "GEXP")
  suffs <- list(nomir)


  sample_data_file <- paste0(dataDir, "/combined_matrix.tsv")
  sample_data <- get_sample_data(sample_data_file)

  labels <- sample_data[[1]]
  dat_raw <- sample_data[[2]]

  message(paste0("loaded sample data\tsamples:", dim(dat_raw)[1], "\tfeatures:",
    dim(dat_raw)[2]))

  message(paste0("num labels: ", length(labels)))

  classification_type <- detect_classification_type(labels)
  message(paste0("detected classification type: ", classification_type))

  message("meth mapper")

  saved_meth_mapper_file <- "meth_mapper_hm450.RData"
  found_meth_mapper_file <- TRUE %in% (list.files(path = homeDir) == saved_meth_mapper_file)

  saved_meth_mapper_file <- paste0(homeDir, "/", saved_meth_mapper_file)
  if (found_meth_mapper_file) {
    message(paste0("found ", saved_meth_mapper_file))
  } else {
    message(paste0("getting a new meth.mapper to save to "), saved_meth_mapper_file)
    hm450 <- get450k()
    meth.mapper <- getNearestGene(hm450)
    save(meth.mapper, file = saved_meth_mapper_file)
  }

  message(paste0("loading meth.mapper from "), saved_meth_mapper_file)
  load(saved_meth_mapper_file)

  meth_map_sample_data(dat_raw, meth.mapper)


  message("load cross-cohort cv feature whitelists")
  cv_feature_whitelists_file <- paste0(dataDir, "/union_of_cv_feature_sets_50.tsv")
  cv_feature_whitelists <- readSetList(cv_feature_whitelists_file)

  message(paste0("num CV feature whitelists: ", length(cv_feature_whitelists)))
  message(names(cv_feature_whitelists))

  tasks <- names(cv_feature_whitelists)
  message("tasks:")
  message(tasks)

  # aklimate expects a dataframe, not a matrix
  dat_df <- as.data.frame(dat_raw)
}


#######################################################

num_of_partitions_per_fold <- 5
num_repeats <- length(tasks)/num_of_partitions_per_fold

if (TRAIN_TEST) {

  repeats_list <- lapply(1:num_repeats, function(x) tasks[seq(num_of_partitions_per_fold *
    (x - 1) + 1, num_of_partitions_per_fold * x, 1)])
  message("repeats_list")
  message(repeats_list)

  message("building reduced models and predicting")
  bacc_results <- foreach(cv_folds = iter(repeats_list), .combine = rbind) %dopar%
    {
      message(cv_folds)
      bacc_for_folds <- foreach(cv_fold = iter(cv_folds), .combine = c) %do%
        {

          message(cv_fold)

          cv_feature_whitelist <- restrict_sample_data_to_feature_whitelist(dat_df,
          cv_feature_whitelists, cv_fold)

          dat_fold <- mlr::createDummyFeatures(dat_df)

          samples_train <- rownames(splits)[splits[, cv_fold] == 0]
          message(paste0("num samples_train: ", length(samples_train)))

          dat_train <- cbind(data.frame(labels = labels[samples_train]),
          dat_fold[samples_train, cv_feature_whitelist, drop = FALSE])

          cutoff_adjusted <- length(cv_feature_whitelist)
          message(paste0("cutoff_adj: ", cutoff_adjusted))

          message("train model")

          # setting importance='impurity_corrected' will compute Gini Importance, or Mean
          # Decrease in Impurity (MDI)
          # https://alexisperrier.com/datascience/2015/08/27/feature-importance-random-forests-gini-accuracy.html
          # vlado recommends importance='permutation', Mean Decrease in Accuracy (MDA)

          # importance_option <- 'none' importance_option <- 'impurity_corrected'
          importance_option <- "permutation"
          message(paste0("using feature importance: ", importance_option))

          rf <- ranger(data = dat_train, dependent.variable.name = "labels",
          always.split.variables = NULL, classification = TRUE, sample.fraction = 0.5,
          num.trees = 3000, mtry = ceiling(cutoff_adjusted/5), min.node.size = 1,
          case.weights = NULL, num.threads = 3, probability = TRUE, respect.unordered.factors = FALSE,
          importance = importance_option, write.forest = TRUE, keep.inbag = TRUE,
          replace = FALSE)

          save(rf, file = paste0(modelsDir, "/", cv_fold, "_", cutoff_adjusted,
          "_cross_cohort_features_rf_reduced_model.RData"))

          if (importance_option != "none") {
          message("save feature importance")
          feature_importance <- rf[["variable.importance"]]
          feature_importance <- feature_importance[order(unlist(feature_importance),
            decreasing = TRUE)]
          feature_importance_df <- as.data.frame(feature_importance)
          feature_importance_file <- paste0(modelsDir, "/", cv_fold, "_",
            cutoff_adjusted, "_cross_cohort_features_rf_reduced_model_feature_importance.tsv")
          write.df(feature_importance_df, row.names.id = "feature", feature_importance_file)
          }

          samples_test <- setdiff(rownames(splits), samples_train)
          message(paste0("num samples_test: ", length(samples_test)))

          dat_test <- dat_df[samples_test, cv_feature_whitelist]
          message(paste0("dim(dat_test):", dim(dat_test)))

          message("predict test samples")
          rf.preds <- predict(rf, dat_test)$predictions

          rownames(rf.preds) <- samples_test

          save(rf.preds, file = paste0(modelsDir, "/", cv_fold, "_", cutoff_adjusted,
          "_cross_cohort_features_rf_reduced_model_predictions.RData"))



          message("compute performance")
          max_labels <- apply(rf.preds[samples_test, ], 1, which.max)
          confusion_matrix_data <- factor(sapply(max_labels, function(x) levels(labels)[x]),
          levels = levels(labels))

          true_labels <- labels[samples_test]

          confM <- caret::confusionMatrix(confusion_matrix_data, true_labels)

          save(confM, file = paste0(modelsDir, "/", cv_fold, "_", cutoff_adjusted,
          "_cross_cohort_features_rf_reduced_model_stats.RData"))



          if (classification_type == "binary") {
          # for binary classification, use this
          bacc_results <- unname(confM$byClass["Balanced Accuracy"])
          } else if (classification_type == "multiclass") {
          # for multiclass classification, use this
          bacc_results <- mean(unname(confM$byClass[, "Balanced Accuracy"]))
          } else {
          message(paste0("**ERROR** CLASSIFICATION_TYPE must be binary or multiclass. CLASSIFICATION_TYPE=",
            classification_type))
          stopifnot(FALSE)
          }

          return(bacc_results)
        }
      mean_bacc_for_folds <- mean(bacc_for_folds)
      message(paste0("mean_bacc_for_folds: ", mean_bacc_for_folds))
      return(bacc_for_folds)

    }
}


#######################################################


message("writing bacc results")
bacc_colnames <- lapply(seq(1, num_of_partitions_per_fold), function(x) paste0("F",
  x))
colnames(bacc_results) <- bacc_colnames

bacc_rownames <- lapply(seq(1, num_repeats), function(x) paste0("R", x))
rownames(bacc_results) <- bacc_rownames

write.df(data.frame(bacc_results, check.names = FALSE), "repeat_fold", paste0(modelsDir,
  "/bacc_cross_cohort_cv_feature_sets.tab"))

message("DONE!")
