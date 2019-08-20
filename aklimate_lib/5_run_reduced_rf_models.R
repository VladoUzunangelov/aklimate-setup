#!/usr/bin/env Rscript

## copyright (c) 2019 Vlado Uzunangelov

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

tasks <- 1:25
message("tasks: ", tasks)

############################################################################

message("load pathways")

# p1 <- readSetList(paste0(featureSetsDir, '/pathcomm_pathways_cleaned')) p1 <-
# p1[!grepl('[[:digit:]]_SMPDB$|Z_SMPDB$', names(p1))]
p1 <- readSetList(paste0(featureSetsDir, "/pathcomm_pathways_cleaned_non_redundant_names.tsv"))
p2 <- readSetList(paste0(featureSetsDir, "/genomic_position_sets.listt"))
p3 <- readSetList(paste0(featureSetsDir, "/genesigdb_human.tab"))
p4 <- readSetList(paste0(featureSetsDir, "/msigdb_c2_c5_no_c2_cp.tab"))

pathways <- c(p1, p2, p3, p4)
max_size_of_pathways <- 1000
pathways <- pathways[sapply(pathways, length) < max_size_of_pathways]


message("sanitize pathway names")

## probably a good idea to sanitize the name a bit
names(pathways) <- gsub("[^\\w\\s]", "_", names(pathways), perl = TRUE)

message(paste("number of pathways to use", length(names(pathways))))

message("load labels")

labels <- as.matrix(read.delim(paste0(cohortDir, "/labels.tsv"), check.names = F,
  stringsAsFactors = F, header = F, row.names = 1))
labels <- factor(labels[, 1])
lbls <- labels


message("load CV fold splits")

splits <- as.matrix(read.delim(paste0(dataDir, "/cv_folds.tsv"), check.names = F,
  stringsAsFactors = F, header = TRUE, row.names = 1))
splits <- splits[, -1]


message("set tasks and seeds")

tasks <- colnames(splits)
tasks <- tasks[1:25]

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

message("building reduced models and predicting")
acc.reduced <- foreach(i = iter(reps.list), .combine = rbind) %dopar% {
  foreach(k = iter(cutoffs), .combine = c) %do% {
    preds <- foreach(j = iter(i), .combine = c) %do% {
      idx.train <- rownames(splits)[splits[, j] == 0]
      idx.test <- setdiff(rownames(splits), idx.train)

      imps <- read.delim(paste0(modelsDir, "/", j, "_aklimate_multiclass_feature_importance.tab"),
        header = TRUE, row.names = 1)
      dat <- mlr::createDummyFeatures(dat)
      k.adj <- min(k, dim(imps)[1])  # just in case fewer important features than required by the cutoff
      rf <- ranger(data = cbind(data.frame(labels = labels[idx.train]), dat[idx.train,
        rownames(imps)[1:k.adj], drop = FALSE]), dependent.variable.name = "labels",
        always.split.variables = NULL, classification = TRUE, sample.fraction = 0.5,
        num.trees = 3000, mtry = ceiling(k.adj/5), min.node.size = 1, case.weights = NULL,
        num.threads = 3, probability = TRUE, respect.unordered.factors = FALSE,
        importance = "none", write.forest = TRUE, keep.inbag = TRUE, replace = FALSE)
      save(rf, file = paste0(modelsDir, "/", j, "_cutoff_", k, "_rf_reduced_model.RData"))
      rf.preds <- predict(rf, dat[idx.test, rownames(imps)[1:k.adj]])$predictions
      rownames(rf.preds) <- idx.test
      save(rf.preds, file = paste0(modelsDir, "/", j, "_cutoff_", k, "_rf_reduced_model_predictions.RData"))

      confM <- caret::confusionMatrix(factor(apply(rf.preds[idx.test, ], 1,
        which.max), levels = levels(labels)), labels[idx.test])
      save(confM, file = paste0(modelsDir, "/", j, "_cutoff_", k, "_rf_reduced_model_stats.RData"))

      mean(unname(confM$byClass[, "Balanced Accuracy"]))
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
