#!/usr/bin/Rscript

# usage, options and doc goes here
argspec <- c("get_aklimate_feature_set_weights.R
Usage:
  Rscript get_aklimate_feature_set_weights.R arg1
Example:
\tRscript get_aklimate_feature_set_weights.R arg1 arg2 arg3 arg4
Options:
\targ1 is a filename to load the AKLIMATE model, such as R1:F1_junkle_final_model.RData.
")

write_feature_set_weights_to_file <- function(feature_set_weights, filename) {
  write.df(data.frame(feature_set_weights), row.names.id = "feature_set_name",
    filename)
}

main <- function(argv) {
  print("argv below:")
  print(argv)
  print("do something")

  message("load full AKLIMATE model")
  load(argv[1])

  message("get model")
  model <- jklm[["junkle.model"]]

  message(paste0("length(model)=", length(model)))

  message("get feature set weights")
  for (i in 1:length(model)) {
    message("i=", i)
    weights_named_vector <- model[[i]][["sorted_kern_weight"]]
    write_feature_set_weights_to_file(weights_named_vector, paste0("AKLIMATE_feature_set_weights_",
      i, ".tsv"))
  }
}

main(commandArgs(TRUE))
