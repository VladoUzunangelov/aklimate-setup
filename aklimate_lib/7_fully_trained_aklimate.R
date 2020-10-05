#!/usr/bin/env Rscript
# vuzunangelov & chrisw 20190606

LOAD_PATHWAYS = TRUE
LOW_EXPRESSION_FILTER = FALSE
QUANTIZE_NUMERIC_DATA = FALSE
LOAD_SAMPLE_DATA_MATRIX = TRUE
LOAD_CV_SPLITS = FALSE

ncpus <- NUMBER_OF_CPUS_TO_USE
stopifnot(is.numeric(ncpus), ncpus > 0, ncpus%%1 == 0)
message("using ", ncpus, " CPUs")

cohortDir <- getwd()
featureSetsDir <- paste0(cohortDir, "/p_store_files")
modelsDir <- paste0(cohortDir, "/models")
dataDir <- paste0(cohortDir, "/data")

homeDir <- featureSetsDir

############################################################################


if (LOAD_PATHWAYS) {

  message("load pathways")

  # using one file of collected pathways.  The names of these pathways have been
  # mapped to simple names.  This was done to avoid the continued problem of name
  # collisions and the like.  collected_pathways_name_mapped.tsv
  p_collected <- readSetList(paste0(featureSetsDir, "/collected_pathways_name_mapped.tsv"))
  pathways <- c(p_collected)

  max_size_of_pathways <- 1000
  pathways <- pathways[sapply(pathways, length) < max_size_of_pathways]

  # need to load the name mappings in order to map back to original pathway names
  p_name_mapping <- read.csv(file = paste0(featureSetsDir, "/pathway_name_mapping.tsv"),
    header = FALSE, sep = "\t", col.names = c("p_id", "p_name"))

  message("sanitize pathway names")

  ## probably a good idea to sanitize the name a bit
  names(pathways) <- gsub("[^\\w\\s]", "_", names(pathways), perl = TRUE)

  message(paste("number of pathways to use", length(names(pathways))))

} else {
  message("skip loading pathways")
}



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

if (LOAD_CV_SPLITS) {

  message("load CV fold splits")

  splits <- as.matrix(read.delim(paste0(dataDir, "/cv_folds.tsv"), check.names = F,
    stringsAsFactors = F, header = TRUE, row.names = 1))
  splits <- splits[, -1]

  message("set tasks and seeds")

  tasks <- colnames(splits)
  # tasks <- tasks[1:25] tasks <- c(tasks[1])

  message("tasks")
  message(tasks)

  # you can change the seeds - this is just to match initial runs seeds <- 11 *
  # (1:length(tasks)) names(seeds) <- tasks

} else {
  message("skip loading cv splits")
  seed <- 11
}


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



  saved_meth_mapper_file <- "meth_mapper_hm450.RData"
  found_meth_mapper_file <- TRUE %in% (list.files(path=homeDir) == saved_meth_mapper_file)

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


message("setup full AKLIMATE training")

set.seed(seed, kind = "L'Ecuyer-CMRG")

classification_type <- CLASSIFICATION_TYPE
message(paste0(" classification_type: ", classification_type))


message("set params for junkle model")


# dat, dat.grp, lbls, fsets, always.add=NULL, rf.pars=list(), junkle.pars=list(),
# store.kernels=FALSE, verbose=FALSE
idx.train <- rownames(dat)
junkle_training_data <- dat[idx.train, ]
junkle_datatypes <- suffs
junkle_training_labels <- labels[idx.train]
junkle_feature_sets <- pathways
junkle_always_add <- NULL

df_for_oob.cv <- data.frame(min.node.prop = 0.01, mtry.prop = 0.25, ntree = 500)

# I got this error:
# [1] "Finished kernel construction"
#Error in 1:which(cumsum(ll) > kcv$pars[kcv$best.id, "nkern"])[1] :
#  NA/NaN argument
#>
#
# To fix it, adjust num_of_trees (not oob.cv) by increasing from the dafault (1000).

num_of_trees = 1000
# num_of_trees = 1500

rf_params <- list(ttype = classification_type, bin.perf = c("bacc"), importance = "permutation",
  min.nfeat = 15, ntree = num_of_trees, sample.frac = 0.5, replace = FALSE, weights = NULL,
  oob.cv = df_for_oob.cv)


# Set nfold based on the size of the smallest class.  We encountered situations
# where the nfold was larger than the size of smallest class.  Best solution is
# to remove the tiny class, but here we attempt to keep all classes.
default_nfold <- 5
num_folds <- default_nfold
for (class in levels(labels)) {
  training_labels <- labels[idx.train]
  training_set_for_class <- (training_labels)[training_labels[] == class]
  num_folds <- min(length(training_set_for_class), num_folds)
}
message("setting num_folds = ", num_folds)

# common warning message:
# The duality gap has been closing very slowly indicating slow convergence.You should examine your kernels for multicollinearity and or change regularization parameters.Alternatively you can increase minIter or decrease tolMinIter.

# MKL solver is converging very slowly.
# address this by using different range of regularization parameter tuning as below.
# The values are log2 scale, so negative values means very close to zero, i.e. low regularization

# small number underflow with rcpp might cause error like this:
# Error in { : task 8 failed - "NA/NaN argument"
# and also many slow convergence messages
# changing lamb range may correct this
junkle_lamb = c(-5, 0)
junkle_lamb = c(-15, 0)

# this is the setting used for TMP
junkle_lamb = c(-8, 3)

junkle_params <- list(topn = 5, nfold = num_folds, subsetCV = TRUE, lamb = junkle_lamb, cvlen = 200, type = "probability")

junkle_store_kernels <- FALSE
junkle_verbose <- TRUE




message("train model")

model_file_name = "junkle_fully_trained_model.RData"
model_file_path = paste0(modelsDir, "/", model_file_name)

if (TRUE %in% (list.files(path=modelsDir) == model_file_name)) {
  message(paste0("loading model from file: ", model_file_path))
  load(model_file_path)
} else {
  jklm <- junkle(junkle_training_data, junkle_datatypes, junkle_training_labels, junkle_feature_sets, junkle_always_add, rf_params, junkle_params, junkle_store_kernels, junkle_verbose)

  save(jklm, file = model_file_path)
}


#######################################################

message("get feature importance")

imps <- rank.features.jklm(jklm)
if (is.null(imps)) {
  message("imps is null!!!")
}
write.df(data.frame(importance = imps), "features", paste0(modelsDir, "/junkle_fully_trained_model_feature_importance.tab"))


#######################################################

message("get feature set weights")


write_feature_set_weights_to_file <- function(feature_set_weights, filename) {
  write.df(data.frame(feature_set_weights), row.names.id = "feature_set_name",
    filename)
}

model <- jklm[["junkle.model"]]
outfilename <- paste0(modelsDir, "/junkle_fully_trained_model_feature_set_importance")

if (classification_type == "binary") {

  weights_named_vector <- model[["sorted_kern_weight"]]
  write_feature_set_weights_to_file(weights_named_vector, outfilename)

} else {

  message(paste0("length(model)=", length(model)))

  for (i in 1:length(model)) {
    message("i=", i)
    weights_named_vector <- model[[i]][["sorted_kern_weight"]]
    write_feature_set_weights_to_file(weights_named_vector, paste0(outfilename,
      "_", i))
  }
}


#######################################################

message("get training predictions")

training_predictions_outfilename <- paste0(modelsDir, "/junkle_fully_trained_model_training_predictions.tsv")
write.df(jklm$preds.train, row.names.id = "sampleID", training_predictions_outfilename)


#######################################################

message("get overall training performance")


jklm.preds <- apply(jklm$preds, 1, function(x) colnames(jklm$preds)[which.max(x)])

confM <- caret::confusionMatrix(factor(jklm.preds, levels = levels(labels)), labels[idx.train])

# get BACC results
if (classification_type == "binary") {
  # for binary classification, use this
  bacc <- unname(confM$byClass["Balanced Accuracy"])
} else if (classification_type == "multiclass") {
  # for multiclass classification, use this
  bacc <- mean(unname(confM$byClass[, "Balanced Accuracy"]))
} else {
  message(paste0(i, " **ERROR** CLASSIFICATION_TYPE must be binary or multiclass. CLASSIFICATION_TYPE=",
    CLASSIFICATION_TYPE))
  stopifnot(FALSE)
}

save(confM, bacc, file = paste0(modelsDir, "/junkle_fully_trained_model_stats_preds.RData"))

#######################################################

message("get subtype training performance")


if (classification_type == "binary") {
  # for binary classification, use this
  subtype_bacc <- confM[["byClass"]]["Balanced Accuracy"]

} else if (classification_type == "multiclass") {
  # for multiclass classification, use this
  subtype_bacc <- confM[["byClass"]][, "Balanced Accuracy"]

} else {
  message(paste0(i, " **ERROR** CLASSIFICATION_TYPE must be binary or multiclass. CLASSIFICATION_TYPE=",
    CLASSIFICATION_TYPE))
  stopifnot(FALSE)
}

subtype_bacc_outfilename <- paste0(modelsDir, "/junkle_fully_trained_model_subtype_bacc.tsv")
write.df(data.frame(subtype_bacc), row.names.id = "class", subtype_bacc_outfilename)


#######################################################

message("DONE !!!")
