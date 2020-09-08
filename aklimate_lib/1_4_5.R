NUMBER_OF_CPUS_TO_USE <- yyy

splits <- as.matrix(read.delim("./data/cv_folds.tsv", check.names = F, stringsAsFactors = F, header = TRUE, row.names = 1))
splits <- splits[, -1]

source("1_load_aklimate_libs.R")
source("4_get_feature_importance_from_full_aklimate.R")
source("5_run_reduced_rf_models.R")
