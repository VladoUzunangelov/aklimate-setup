# A toy example for running AKLIMATE on a binary classification task

# vuzunangelov & chrisw 20210210

message("==> load AKLIMATE library")

library(aklimate)

message("==> finish load libraries")

############################################

message("==> begin load data")

message("==> setup dirs")
workDir <- getwd()

task.dir <- paste(workDir, "example_results", sep = "/")
if (!dir.exists(task.dir)) dir.create(task.dir)

# tasks <- c(1)

rdata_filename <- "binary_classification_example_data.rdata"
example_data_rdata_file <- paste0(workDir, "/", rdata_filename)

if (TRUE %in% (list.files(path = workDir) == "binary_classification_example_data.rdata")) {
  message(paste0("==> load data from rdata file: ", example_data_rdata_file))
  load(example_data_rdata_file)

  pathways <- binary_classification_example_data$pathways
  labels <- binary_classification_example_data$labels
  sample_data <- binary_classification_example_data$sample_data
  train_test_split <- binary_classification_example_data$train_test_split

} else {
  message("==> load data from individual data files")

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

  num_cv_folds <- 4
  tasks <- 1:num_cv_folds

  splits <- createFolds(labels[, 1], k = num_cv_folds, list = TRUE, returnTrain = TRUE)

  train_test_split <- splits[[1]]



  message("==> finish load data")


  message("==> generate example data object")
  binary_classification_example_data <- list(pathways = pathways, labels = labels,
    sample_data = sample_data, splits = splits, tasks = tasks, train_test_split = train_test_split)

  class(binary_classification_example_data) <- "aklimate_example_data"

  save(binary_classification_example_data, file = example_data_rdata_file)

}

#######################

message("==> begin running AKLIMATE")

set.seed(11, kind = "L'Ecuyer-CMRG")

message("==> setting parameters. Some settings used here sacrifice classification performance for speed.")
message("description of parameters at https://github.com/VladoUzunangelov/aklimate")

idx.train <- rownames(labels)[train_test_split]
idx.test <- setdiff(rownames(labels), idx.train)

training_data <- sample_data[idx.train, ]
datatypes <- list(exp = c("exp"))
training_lables <- labels[idx.train, 1]
feature_sets <- pathways
global_features_names <- NULL

rf_params <- list(ttype = "binary", bin.perf = "bacc", regression.q = 0.05, importance = "impurity",
  min.nfeat = 15, ntree = 2000, sample.frac = 0.5, replace = FALSE, weights = NULL,
  oob.cv = data.frame(min.node.prop = 0.01, mtry.prop = 0.25, ntree = 500))

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

message("==> test model")
tc <- Sys.time()
aklimate_model.preds <- predict(aklimate_model, sample_data, pathways, NULL, FALSE)$preds

td <- Sys.time()

t_test <- difftime(td, tc, units = "secs")
message(paste0("testing time: ", t_test))

message(paste0("train + test time: ", t_train + t_test))

message("==> get ranked features")
ranked_features <- rank_features(akl_obj=aklimate_model)

message("==> get relative importance of datatypes")
datatype_importance <- rank_importance_type(suffs=c("exp"), ranked_features)


message("==> finish running AKLIMATE")
