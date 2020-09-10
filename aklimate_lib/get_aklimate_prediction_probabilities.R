#!/usr/bin/Rscript

# usage, options and doc goes here
argspec <- c("get_aklimate_prediction_probabilities.R
Usage:
    Rscript get_aklimate_prediction_probabilities.r arg1 arg2 arg3 arg4
Example:
\tRscript get_aklimate_prediction_probabilities.r arg1 arg2 arg3 arg4
Options:
\targ1 must be full or reduced to indicate full AKLIMATE model or reduced AKLIMATE model.
\targ2 must be test or train to indicate test samples or training samples.
\targ3 is a filename. as follows:
\t\tfull test - XXX_junkle_final_model_stats_preds.RData
\t\tfull train - XXX_junkle_final_model.RData
\t\treduced test - XXX_reduced_model_predictions.RData
\t\treduced train - XXX_reduced_model.RData
\targ4 is a filename required just for retrieving sampleIDs for reduced model training set probs:
\t\treduced train - XXX_junkle_final_model.RData OR XXX_idx_train_labels.RData
")


get_full_model_test_probs <- function(filename, outFile= "full_test_probs.tsv") {
    load(filename)
    #write.df(data.frame(as.list(jklm.preds)),row.names.id = "sampleID", outFile)
    write.df(data.frame(as.list(jklm.preds)),row.names.id = "sampleID", outFile)
}


get_full_model_train_probs <- function(filename, outFile = "full_train_probs.tsv") {
    load(filename)
    write.df(jklm$preds.train, row.names.id = "sampleID", outFile)
}


get_reduced_model_test_probs <- function(filename, outFile = "reduced_test_probs.tsv") {
    load(filename)
    write.df(rf.preds,row.names.id = "sampleID", outFile)
}

get_reduced_model_train_probs <- function(filename1, filename2, outFile = "reduced_train_probs.tsv") {
  load(filename1)
  load(filename2)

  if (exists("training_sample_labels")) {

    # for debugging.... checking if prediction matches true label z <-
    # sapply(training_sample_labels, function(x) { levels(training_sample_labels)[x]
    # }) a <- cbind(rf$predictions, z)

    ids <- names(training_sample_labels)

  } else {
    ids <- jklm$idx.train
  }
  pred <- data.frame(rf$predictions)
  rownames(pred) <- ids

  write.df(pred, row.names.id = "sampleID", outFile)
}


main <- function(argv) {
    print("argv below:")
    print(argv)
    print("do something")
    if (argv[1] == "full" && argv[2] == "test") {
      print("get full model test probs")
      get_full_model_test_probs(argv[3])
    } else if (argv[1] == "full" && argv[2] == "train") {
      print("get full model train probs")
      get_full_model_train_probs(argv[3])
    } else if (argv[1] == "reduced" && argv[2] == "test") {
      print("get reduced model test probs")
      get_reduced_model_test_probs(argv[3])
    } else if (argv[1] == "reduced" && argv[2] == "train") {
      print("get reduced model train probs")
      if(length(argv) < 4){
        print("ERROR: 2 files required for reduced train")
        stopifnot(FALSE)
      }
      get_reduced_model_train_probs(argv[3], argv[4])
    } else {
      print("ERROR: arg1 must be full or reduced. arg2 must be test or train.")
    }
}

main(commandArgs(TRUE))
