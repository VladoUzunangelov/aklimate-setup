#!/usr/bin/Rscript
# 20220407 chris
# 20210708 chris

# usage, options and doc goes here
argspec <- c("10_predict_samples_using_reduced_models.R
Usage:
  Rscript 10_predict_samples_using_reduced_models.R arg1 arg2
Example:
\tRscript 10_predict_samples_using_reduced_models.R arg1 arg2
Options:
\targ1 cohort (ACC, BLCA, BRCA, CESC, COADREAD, ESCC, GEA, HNSC, KIRCKICH, KIRP, LGGGBM, LIHCCHOL, LUAD, LUSC, MESO, OV, PAAD, PCPG, PRAD, SARC, SKCM, TGCT, THCA, THYM, UCEC, UVM)
\targ2 datatype (GEXP, CNVR, METH, MULTI, TOP)
\targ3 sample_data (Tab-separated matrix file. Each row is a sample. Each column is a feature.)")

# The dev environment used:
# R version 3.4.4 (2018-03-15)
# ranger_0.11.2
# doParallel_1.0.14

library(doParallel)
registerDoParallel(cores = 1)

library(ranger)


main <- function(argv) {
  cohort <- argv[1]
  datatype <- argv[2]
  sample_data_file <- argv[3]

  message(paste0("cohort: ", cohort))
  message(paste0("datatype: ", datatype))
  message(paste0("sample_data_file: ", sample_data_file))

  cohort_datatype <- toupper(paste0(cohort, "_", datatype))

  sample_data_dir <- dirname(sample_data_file)
  predictions_file_name <- paste0(cohort_datatype, "_predictions_for_", basename(sample_data_file))
  predictions_file_path <- paste0(sample_data_dir, "/", predictions_file_name)

  # data_dir is where R will find model RData files.
  data_dir <- "/aklimate"
  # sample_data_file <- paste0(data_dir, '/', cohort, '_combined_matrix.tsv')

  # file_search_string <- paste0(cohort, "_*_junkle_fully_trained_reduced_model.RData")
  file_search_string <- paste0(cohort_datatype, "_rf_reduced_model.RData")
  glob_files <- Sys.glob(file.path(data_dir, file_search_string))

  num_results <- length(glob_files)

  if (num_results == 0) {
    message(paste0("ERROR: did not find any model file for ", cohort_datatype))
  } else if (num_results == 1) {
    message(paste0("found a model file for predicting subtypes for ", cohort_datatype))
    model_file_path <- glob_files[1]
  } else if (num_results > 1) {
    message(paste0("ERROR: found multiple model files for ", cohort_datatype))
  }

  stopifnot(exists("model_file_path"))

  # message(paste0('model_file_path: ', model_file_path))

  # load rf
  load(model_file_path)

  model_features <- names(rf$variable.importance)
  message(paste0("The model uses ", length(model_features), " features:"))
  message(paste0(" ", model_features))

  saved_meth_mapper_file <- paste0(data_dir, "/", "meth_mapper_hm450_20200416.RData")

  #############################

  message(paste0("load sample data from: ", sample_data_file))

  nomir = c("MUTA", "CNVR", "METH", "GEXP")
  suffs <- list(nomir)

  dat <- read.delim(sample_data_file, header = T, row.names = 1, check.names = FALSE)

  message("  feature matrix dimensions:")
  message(paste0(" ", dim(dat)))


  message("processing sample data")

  # load meth.mapper
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


  #######################################################

  message("predict subtypes")

  pred_results <- predict(rf, dat)

  predictions <- pred_results$predictions
  rownames(predictions) <- rownames(dat)

  max_labels <- apply(predictions, 1, which.max)

  labels <- colnames(predictions)

  max_labels <- sapply(max_labels, function(x) labels[x])
  names(max_labels) <- rownames(predictions)


  predictions_final_labels_and_scores <- cbind(max_labels, predictions)

  write.table(predictions_final_labels_and_scores, file = predictions_file_path,
    row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)

  message(paste0("Done! Predictions saved to: ", predictions_file_path))
}


main(commandArgs(TRUE))
