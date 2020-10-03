#!/usr/bin/env Rscript

## copyright (c) 2019 Vlado Uzunangelov

# NUMBER_OF_CPUS_TO_USE = 28
GENERATE_FEATURE_IMPORTANCE_FILES = TRUE

message("Set some variables, load some libraries, source some paths.")

ncpus <- NUMBER_OF_CPUS_TO_USE
stopifnot(is.numeric(ncpus), ncpus > 0, ncpus%%1 == 0)
message("using ", ncpus, " CPUs")

cohortDir <- getwd()
modelsDir <- paste0(cohortDir, "/models")

# tasks <- 1:25
tasks <- colnames(splits)[1:25]
message("tasks: ", tasks)


#######################################################

if (GENERATE_FEATURE_IMPORTANCE_FILES) {
  message("getting feature importance from full models")
  z1 <- foreach(i = iter(tasks)) %dopar% {

    full_model_importance_file = paste0(i, "_aklimate_multiclass_feature_importance.tab")
    full_model_importance_path = paste0(modelsDir, "/", full_model_importance_file)

    if (TRUE %in% (list.files(path=modelsDir) == full_model_importance_file)) {
      message(paste0("found feature importance file: ", full_model_importance_file))
    } else {
      message(paste0("getting feature importance for: ", i))
      load(paste0(modelsDir, "/", i, "_junkle_final_model.RData"))
      if (is.null(jklm)) {
        message("jklm object is null for ", i)
      }
      imps <- rank.features.jklm(jklm)
      if (is.null(imps)) {
        messge("imps is null for ", i)
      }
      write.df(data.frame(importance = imps), "features", full_model_importance_path)
    }
  }
} else {
  message("skipping generating importance files")
}

#######################################################

message("DONE!")
