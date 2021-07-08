#!/usr/bin/Rscript
# 20210708 chrisw

# usage, options and doc goes here
argspec <- c("script.r
Usage:
  Rscript script.r arg1 arg2 arg3
Example:
\tRscript script.r arg1 arg2 arg3
Options:
\targ1 cohort
\targ2 datatype
\targ3 feature_list_size")

NUMBER_OF_CPUS_TO_USE <- 9
ncpus <- NUMBER_OF_CPUS_TO_USE
stopifnot(is.numeric(ncpus), ncpus > 0, ncpus%%1 == 0)
message("using ", ncpus, " CPUs")

library(doParallel)
registerDoParallel(cores = ncpus)

library(ranger)
library(FDb.InfiniumMethylation.hg19)

LOAD_PATHWAYS = FALSE
LOW_EXPRESSION_FILTER = FALSE
QUANTIZE_NUMERIC_DATA = FALSE
LOAD_SAMPLE_DATA_MATRIX = TRUE
LOAD_CV_SPLITS = FALSE

main <- function(argv) {
  message("argv below:")
  message(argv)
  # print('do something')

  cohort <- argv[1]
  datatype <- argv[2]
  feature_list_size <- as.numeric(argv[3])

  data_dir <- "./reduced_model_input"
  sample_data_file <- paste0(data_dir, "/", cohort, "_combined_matrix.tsv")
  feature_list_file <- paste0(data_dir, "/", cohort, "_", datatype, "_junkle_fully_trained_model_feature_importance.tab")

  output_dir <- "./models"
  model_output_file <- paste0(cohort, "_", datatype, "_", feature_list_size, "_junkle_fully_trained_reduced_model.RData")
  reduced_model_feature_importance_file <- paste0(cohort, "_", datatype, "_", feature_list_size,
    "_rf_reduced_model_feature_importance.tsv")

  saved_meth_mapper_file <- "meth_mapper_hm450_20200416.RData"

  #############################

  if (LOAD_SAMPLE_DATA_MATRIX) {
    message(paste0("load sample data from: ", sample_data_file))

    nomir = c("MUTA", "CNVR", "METH", "GEXP")
    suffs <- list(nomir)

    dat <- read.delim(sample_data_file, header = T, row.names = 1, check.names = FALSE)

    message("  combined matrix dimensions:")
    message(paste0(" ", dim(dat)))

    message("extract labels")

    # labels in column 2
    labels <- factor(dat[, 1])
    names(labels) <- rownames(dat)

    num_labels <- length(levels(labels))
    if (num_labels == 2) {
      CLASSIFICATION_TYPE = "binary"
    } else if (num_labels > 2) {
      CLASSIFICATION_TYPE = "multiclass"
    } else {
      message(paste0("ERROR: detected ", num_labels, " labels"))
      stopifnot(FALSE)
    }

    message(CLASSIFICATION_TYPE)
    message(paste0("number of samples:", length(labels)))
    message("  labels:")
    message(paste0(" ", levels(labels)))


    # remove column 2, the labels
    dat <- dat[, -1]

    message("  feature matrix dimensions:")
    message(paste0(" ", dim(dat)))


    message("meth mapper")


    found_meth_mapper_file <- TRUE %in% (list.files(path = data_dir) == saved_meth_mapper_file)

    saved_meth_mapper_file <- paste0(data_dir, "/", saved_meth_mapper_file)
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

    updated.names <- foreach(i = iter(colnames(dat)), .combine = c) %dopar% {
      if (grepl("METH", i, ignore.case = TRUE)) {
        cg <- strsplit(i, ":")[[1]][4]
        gsub(paste0(cg, "::"), paste0(cg, ":", meth.mapper[cg, 4], ":"),
          i)
      } else {
        i
      }
    }

    colnames(dat) <- updated.names

    if (LOW_EXPRESSION_FILTER) {
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
      exp_features_upper75_by_sd = names(exp_sds)[exp_sds > quantile(exp_sds,
        probs = 0.25)]
      exp_all_features <- colnames(exp_matrix)
      exp_keep_features <- intersect(exp_features_upper75_by_mean, exp_features_upper75_by_sd)
      exp_drop_features <- setdiff(exp_all_features, exp_keep_features)

      full_keep_features <- setdiff(colnames(dat), exp_drop_features)

      message(paste("number of expression features to drop", length(exp_drop_features)))

      exp_matrix <- NULL

      message(paste("number of features in full data matrix", length(colnames(dat))))

      dat <- dat[sampleIDs, full_keep_features]
      message(paste("number of features in new data matrix", length(colnames(dat))))

    } else {
      message("skip low expression filter")
    }

    if (QUANTIZE_NUMERIC_DATA) {
      message("quantize numeric datatypes")

      dat <- foreach(datatype_str = iter(nomir[!grepl("MUTA|CNVR", nomir)]),
        .combine = cbind) %do% {
        message(paste("quantize ", datatype_str))
        quantize.data(dat[sampleIDs, grepl(paste0("N:", datatype_str, ":"),
          colnames(dat))], nbr = 5, idx = sampleIDs)
      }
    } else {
      message("skip quantizing numeric datatypes")
    }

    # aklimate expects a dataframe, not a matrix
    dat <- as.data.frame(dat)

  } else {
    message("skip loading sample data matrix")
  }

  #######################################################

  message(paste("load feature importance scores from: ", feature_list_file))

  imps <- read.delim((feature_list_file), header = TRUE, row.names = 1)

  # reverse sort importances
  a <- imps[order(-imps$importance), ]
  names(a) <- rownames(imps)[order(-imps$importance)]
  imps <- as.data.frame(a)
  names(imps) <- "importance"

  message("  head of feature importance scores: ")
  message(paste0(" ", head(imps)))

  message("  dim of feature importance scores: ")
  message(paste0(" ", dim(imps)))

  #######################################################

  message("train reduced model")



  model_file_path <- paste0(output_dir, "/", model_output_file)
  if (TRUE %in% (list.files(path = output_dir) == model_output_file)) {
    message(paste0("loading model from file: ", model_file_path))
    load(model_file_path)
  } else {
    message("train a new model")

    seed <- 11
    set.seed(seed, kind = "L'Ecuyer-CMRG")

    dat <- mlr::createDummyFeatures(dat)

    k.adj <- min(feature_list_size, dim(imps)[1])  # just in case fewer important features than required by the cutoff
    message(paste0("k.adj: ", k.adj))

    idx.train <- rownames(dat)
    message(paste0("length(idx.train): ", length(idx.train)))

    training_feature_ids <- rownames(imps)[1:k.adj]

    training_labels <- labels[idx.train]
    training_feature_data <- dat[idx.train, training_feature_ids, drop = FALSE]

    training_data <- cbind(data.frame(labels = labels[idx.train]), training_feature_data)

    message("  head of training_data column names:")
    message(paste0(" ", head(colnames(training_data))))

    rf <- ranger(data = training_data, dependent.variable.name = "labels", always.split.variables = NULL,
      classification = TRUE, sample.fraction = 0.5, num.trees = 3000, mtry = ceiling(k.adj/5),
      min.node.size = 1, case.weights = NULL, num.threads = 3, probability = TRUE,
      respect.unordered.factors = FALSE, importance = "permutation", write.forest = TRUE,
      keep.inbag = TRUE, replace = FALSE)

    save(rf, file = model_file_path)

    message(paste0("trained model written to: ", model_file_path))
  }

  reduced_model_feature_importance_file_path <- paste0(output_dir, "/", reduced_model_feature_importance_file)
  message(paste0("write feature importance from reduced RF model to file", reduced_model_feature_importance_file_path))
  colnames <- paste0("featureID\tfeature_importance_score_for_", cohort, "_", datatype,
    "_", feature_list_size)
  write.table(rf$variable.importance, file = reduced_model_feature_importance_file_path,
    quote = FALSE, sep = "\t", col.names = c(colnames))

  message("DONE !!!")
}


main(commandArgs(TRUE))
