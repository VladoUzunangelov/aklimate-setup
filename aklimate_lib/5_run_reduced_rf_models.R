#!/usr/bin/env Rscript

## copyright (c) 2019 Vlado Uzunangelov

# CLASSIFICATION_TYPE <- "binary"
# CLASSIFICATION_TYPE <- "multiclass"
# NUMBER_OF_CPUS_TO_USE = 28
LOW_EXPRESSION_FILTER = FALSE
QUANTIZE_NUMERIC_DATA = FALSE
LOAD_SAMPLE_DATA_MATRIX = TRUE

message("Set some variables, load some libraries, source some paths.")

ncpus <- NUMBER_OF_CPUS_TO_USE
stopifnot(is.numeric(ncpus), ncpus > 0, ncpus%%1 == 0)
message("using ", ncpus, " CPUs")

cohortDir <- getwd()
featureSetsDir <- paste0(cohortDir, "/p_store_files")
modelsDir <- paste0(cohortDir, "/models")
dataDir <- paste0(cohortDir, "/data")
homeDir <- featureSetsDir

############################################################################


message("load pathways")

# using one file of collected pathways.
# The names of these pathways have been mapped to simple names.
# This was done to avoid the continued problem of name collisions and the like.
# collected_pathways_name_mapped.tsv
p_collected <- readSetList(paste0(featureSetsDir, "/collected_pathways_name_mapped.tsv"))
pathways <- c(p_collected)

max_size_of_pathways <- 1000
pathways <- pathways[sapply(pathways, length) < max_size_of_pathways]

# need to load the name mappings in order to map back to original pathway names
p_name_mapping <- read.csv(file=paste0(featureSetsDir, "/pathway_name_mapping.tsv"), header = FALSE, sep = "\t", col.names = c("p_id", "p_name"))

message("sanitize pathway names")

## probably a good idea to sanitize the name a bit
names(pathways) <- gsub("[^\\w\\s]", "_", names(pathways), perl = TRUE)

message(paste("number of pathways to use", length(names(pathways))))



message("load labels")

labels <- as.matrix(read.delim(paste0(cohortDir, "/labels.tsv"), check.names = F,
  stringsAsFactors = F, header = F, row.names = 1))
labels <- factor(labels[, 1])
lbls <- labels

num_labels <- length(levels(labels))
if (num_labels == 2) {
  CLASSIFICATION_TYPE = "binary"
} else if (num_labels > 2) {
  CLASSIFICATION_TYPE = "multiclass"
} else {
  message(paste0("ERROR: detected ", num_labels, " labels"))
  stopifnot(FALSE)
}

message("load CV fold splits")

splits <- as.matrix(read.delim(paste0(dataDir, "/cv_folds.tsv"), check.names = F,
  stringsAsFactors = F, header = TRUE, row.names = 1))
splits <- splits[, -1]



message("set tasks and seeds")

tasks <- colnames(splits)
tasks <- tasks[1:25]
#tasks <- c(tasks[1])

message("tasks")
message(tasks)

# you can change the seeds - this is just to match initial runs seeds <- 11 *
# (1:length(tasks)) names(seeds) <- tasks

#############################


if (LOAD_SAMPLE_DATA_MATRIX) {
  message("load sample data")

  ## MIR is also there - about 800 MIRs right now probably best to not include them
  ## - I have a scheme to indluce them by putting a MIR in every pathway that has
  ## one of its targets, but it didn't work well when I tested it
  nomir = c("MUTA", "CNVR", "METH", "GEXP")
  suffs <- list(nomir)

  dat <- read.delim(paste0(dataDir, "/combined_matrix.tsv"), header = T, row.names = 1,
    check.names = FALSE)
  dat <- dat[, -1]


  message("meth mapper")

  saved_meth_mapper_file <- "meth_mapper_hm450_20200416.RData"
  found_meth_mapper_file <- TRUE %in% (list.files(path=homeDir) == saved_meth_mapper_file)

  saved_meth_mapper_file <- paste0(homeDir, "/meth_mapper_hm450_20200416.RData")
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
      gsub(paste0(cg, "::"), paste0(cg, ":", meth.mapper[cg, 4], ":"), i)
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

  } else {
    message("skip low expression filter")
  }

  if (QUANTIZE_NUMERIC_DATA) {
    message("quantize numeric datatypes")

    dat <- foreach(datatype_str = iter(nomir[!grepl("MUTA|CNVR", nomir)]), .combine = cbind) %do%
      {
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

cutoffs <- c(5, 10, 20, 50, 100, 200, 500, 1000, 1500)
message("using reduced model sizes: ", cutoffs)

reps.list <- lapply(1:5, function(x) tasks[seq(5 * (x - 1) + 1, 5 * x, 1)])
#reps.list <- reps.list[1]
#reps.list[[1]] <- c("R1:F1")
message("reps.list")
message(reps.list)

message("building reduced models and predicting")
acc.reduced <- foreach(i = iter(reps.list), .combine = rbind) %dopar% {
  foreach(k = iter(cutoffs), .combine = c) %do% {
    preds <- foreach(j = iter(i), .combine = c) %do% {

      stats_file_name = paste0(j, "_cutoff_", k, "_rf_reduced_model_stats.RData")
      preds_file_name = paste0(j, "_cutoff_", k, "_rf_reduced_model_predictions.RData")
      model_file_name = paste0(j, "_cutoff_", k, "_rf_reduced_model.RData")
      feature_importance_file_name = paste0(j, "_aklimate_multiclass_feature_importance.tab")

      stats_file_path = paste0(modelsDir, "/", stats_file_name)
      preds_file_path = paste0(modelsDir, "/", preds_file_name)
      model_file_path = paste0(modelsDir, "/", model_file_name)
      feature_importance_file_path = paste0(modelsDir, "/", feature_importance_file_name)

	  REDUCED_RF_FEATURE_IMPORTANCE_FILE_NAME = paste0(j, "_cutoff_", k, "_rf_reduced_model_feature_importance.tsv")
	  REDUCED_RF_FEATURE_IMPORTANCE_FILE_PATH = paste0(modelsDir, "/", REDUCED_RF_FEATURE_IMPORTANCE_FILE_NAME)

      idx.train <- rownames(splits)[splits[, j] == 0]
      idx.test <- setdiff(rownames(splits), idx.train)

      imps <- read.delim((feature_importance_file_path), header = TRUE, row.names = 1)

      dat <- mlr::createDummyFeatures(dat)

      k.adj <- min(k, dim(imps)[1])  # just in case fewer important features than required by the cutoff

      if (TRUE %in% (list.files(path=modelsDir) == model_file_name)) {
        message(paste0("loading model from file: ", model_file_path))
        load(model_file_path)
      } else {
        rf <- ranger(data = cbind(data.frame(labels = labels[idx.train]), dat[idx.train,
          rownames(imps)[1:k.adj], drop = FALSE]), dependent.variable.name = "labels",
          always.split.variables = NULL, classification = TRUE, sample.fraction = 0.5,
          num.trees = 3000, mtry = ceiling(k.adj/5), min.node.size = 1, case.weights = NULL,
          num.threads = 3, probability = TRUE, respect.unordered.factors = FALSE,
          importance = "permutation", write.forest = TRUE, keep.inbag = TRUE, replace = FALSE)

        save(rf, file = model_file_path)
      }

      if (TRUE %in% (list.files(path=modelsDir) == preds_file_name)) {
        message(paste0("loading preds from file: ", preds_file_path))
        load(preds_file_path)
      } else {
        rf.preds <- predict(rf, dat[idx.test, rownames(imps)[1:k.adj]])$predictions
        rownames(rf.preds) <- idx.test

        save(rf.preds, file = preds_file_path)
      }

      if (TRUE %in% (list.files(path=modelsDir) == stats_file_name)) {
        message(paste0("loading stats from file: ", stats_file_path))
        load(stats_file_path)
      } else {
        # confM <- caret::confusionMatrix(factor(apply(rf.preds[idx.test, ], 1,
        # which.max), levels = levels(labels)), labels[idx.test])

        # cm_data_old <- factor(apply(rf.preds[idx.test, ], 1, which.max), levels =
        # levels(labels))

        max_labels <- apply(rf.preds[idx.test, ], 1, which.max)
        cm_data <- factor(sapply(max_labels, function(x) levels(labels)[x]),
          levels = levels(labels))

        cm_true_labels <- labels[idx.test]

        confM <- caret::confusionMatrix(cm_data, cm_true_labels)

        save(confM, file = stats_file_path)
      }


      classification_type <- CLASSIFICATION_TYPE

      if (classification_type == "binary") {
        # for binary classification, use this
        unname(confM$byClass["Balanced Accuracy"])
      } else if (classification_type == "multiclass") {
        # for multiclass classification, use this
        mean(unname(confM$byClass[, "Balanced Accuracy"]))
      } else {
        message(paste0("**ERROR** CLASSIFICATION_TYPE must be binary or multiclass. CLASSIFICATION_TYPE=",
          CLASSIFICATION_TYPE))
        stopifnot(FALSE)
      }
	  
	  message("write feature importance from reduced RF model to file")
	  colnames<-paste0("featureID	feature_importance_score_for_", j, "_cutoff_", k)
	  write.table(rf$variable.importance, file=REDUCED_RF_FEATURE_IMPORTANCE_FILE_PATH, quote=FALSE, sep="\t", col.names=c(colnames))

      ## mean(unname(confM$overall['Accuracy']))
    }

    mean(preds)

  }
}





message("writing accuracy of reduced feature sets")
colnames(acc.reduced) <- cutoffs
write.df(data.frame(acc.reduced, check.names = FALSE), "repetition", paste0(modelsDir,
  "/balanced_accuracy_reduced_feature_sets.tab"))

######################################################

message("write overall stats for reduced models")
acc.reduced <- foreach(i = iter(reps.list), .combine = rbind) %dopar% {
  foreach(k = iter(cutoffs), .combine = c) %do% {
    preds <- foreach(j = iter(i), .combine = c) %do% {
      load(paste0(modelsDir, "/", j, "_cutoff_", k, "_rf_reduced_model_stats.RData"))
      mean(unname(confM$overall["Accuracy"]))
    }

    mean(preds)

  }
}

colnames(acc.reduced) <- cutoffs
write.df(data.frame(acc.reduced, check.names = FALSE), "repetition", paste0(modelsDir,
  "/accuracy_reduced_feature_sets.tab"))

message("DONE!")
