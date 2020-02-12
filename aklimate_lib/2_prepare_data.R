#!/usr/bin/env Rscript
# vuzunangelov & chrisw 20190606

## Run the commands in this script after those in 1_load_aklimate_libs.R.

QUANTIZE_NUMERIC_DATA = FALSE

message("load sample data")

## MIR is also there - about 800 MIRs right now probably best to not include them
## - I have a scheme to indluce them by putting a MIR in every pathway that has
## one of its targets, but it didn't work well when I tested it
nomir = c("MUTA", "CNVR", "METH", "GEXP")
suffs <- list(nomir)

dat <- read.delim("./data/combined_matrix.tsv", header = T, row.names = 1, check.names = FALSE)
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


message("filter out low expression features")
# Theo has named all gene expression features with 'N:GEXP:'

sampleIDs <- rownames(dat)
exp_matrix <- dat[sampleIDs, grepl(paste0("N:GEXP:"), colnames(dat))]

message(paste("number of features in input expression matrix", length(colnames(exp_matrix))))

exp_means <- colMeans(exp_matrix)
exp_sds <- apply(exp_matrix, 2, sd)

## exclude genes that are in the bottom quartile by mean or sd
exp_features_upper75_by_mean = names(exp_means)[exp_means > quantile(exp_means, probs = 0.25)]
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


if (QUANTIZE_NUMERIC_DATA) {
  message("quantize numeric datatypes")

  dat <- foreach(datatype_str = iter(nomir[!grepl("MUTA|CNVR", nomir)]), .combine = cbind) %do%
    {
      message(paste("quantize ", datatype_str))
      quantize.data(dat[sampleIDs, grepl(paste0("N:", datatype_str, ":"), colnames(dat))],
        nbr = 5, idx = sampleIDs)
    }
}

# aklimate expects a dataframe, not a matrix
dat <- as.data.frame(dat)

message("load pathways")

homeDir <- "./p_store_files"
workDir <- "./models/"

# p1 <- readSetList(paste0(homeDir, '/pathcomm_pathways_cleaned')) p1 <-
# p1[!grepl('[[:digit:]]_SMPDB$|Z_SMPDB$', names(p1))]

# using individual pathway files
# p1 <- readSetList(paste0(homeDir, "/pathcomm_pathways_cleaned_non_redundant_names.tsv"))
# p2 <- readSetList(paste0(homeDir, "/genomic_position_sets.listt"))
# p3 <- readSetList(paste0(homeDir, "/genesigdb_human.tab"))
# p4 <- readSetList(paste0(homeDir, "/msigdb_c2_c5_no_c2_cp.tab"))
#
# pathways <- c(p1, p2, p3, p4)

# using one file of collected pathways.
# The names of these pathways have been mapped to simple names.
# This was done to avoid the continued problem of name collisions and the like.
# collected_pathways_name_mapped.tsv
p_collected <- readSetList(paste0(homeDir, "/collected_pathways_name_mapped.tsv"))
pathways <- c(p_collected)

max_size_of_pathways <- 1000
pathways <- pathways[sapply(pathways, length) < max_size_of_pathways]

# need to load the name mappings in order to map back to original pathway names
p_name_mapping <- read.csv(file=paste0(homeDir, "/pathway_name_mapping.tsv"), header = FALSE, sep = "\t", col.names = c("p_id", "p_name"))

message("sanitize pathway names")

## probably a good idea to sanitize the name a bit
names(pathways) <- gsub("[^\\w\\s]", "_", names(pathways), perl = TRUE)

message(paste("number of pathways to use", length(names(pathways))))

message("load labels")

labels <- as.matrix(read.delim("./labels.tsv", check.names = F, stringsAsFactors = F,
  header = F, row.names = 1))
labels <- factor(labels[, 1])


message("load CV fold splits")

splits <- as.matrix(read.delim("./data/cv_folds.tsv", check.names = F, stringsAsFactors = F,
  header = TRUE, row.names = 1))
splits <- splits[, -1]


message("set tasks and seeds")

tasks <- colnames(splits)

# you can change the seeds - this is just to match initial runs
seeds <- 11 * (1:length(tasks))
names(seeds) <- tasks



message("define worker.f")

worker.f <- function(tasks) {

  res <- foreach(i = iter(tasks)) %do% {

    set.seed(seeds[i], kind = "L'Ecuyer-CMRG")

    idx.train <- rownames(splits)[splits[, i] == 0]
    idx.test <- setdiff(rownames(splits), idx.train)

  # classification_type <- "binary"
  # classification_type <- "multiclass"
    classification_type <- CLASSIFICATION_TYPE
    message(paste0("classification_type: ", classification_type))


    args_list_1 <- list(ttype = classification_type, bin.perf = c("bacc"), importance = "permutation",
      min.nfeat = 15, ntree = 1000, sample.frac = 0.5, replace = FALSE, weights = NULL,
      oob.cv = data.frame(min.node.prop = 0.01, mtry.prop = 0.25, ntree = 500))

    args_list_2 <- list(topn = 5, subsetCV = TRUE, lamb = c(-8, 3), cvlen = 200, type = "probability")

    # train AKLIMATE model
	# works with Spicer.R with md5sum: 71b8fbb55747b623d565960ff7737f73
	# does not work with Spicer.R with md5sum: c1fffbc913a9a82571d6a00116c6dadd
    jklm <- junkle(dat[idx.train, ], suffs, labels[idx.train], pathways, NULL, args_list_1,
      args_list_2, FALSE, TRUE)

    save(jklm, file = paste0(workDir, "/", i, "_junkle_final_model.RData"))

    jklm.preds <- predict.junkle(jklm, dat[c(idx.train, idx.test), ], pathways, NULL,
      FALSE)$preds
    jklm.preds <- apply(jklm.preds, 1, function(x) colnames(jklm.preds)[which.max(x)])

    confM <- caret::confusionMatrix(factor(jklm.preds, levels = levels(labels)), labels[idx.test])

    # get BACC results
    if (classification_type == "binary") {
      # for binary classification, use this
      bacc <- unname(confM$byClass["Balanced Accuracy"])
    } else if (classification_type == "multiclass") {
      # for multiclass classification, use this
      bacc <- mean(unname(confM$byClass[, "Balanced Accuracy"]))
    } else {
      message(paste0("**ERROR** CLASSIFICATION_TYPE must be binary or multiclass. CLASSIFICATION_TYPE=",
        CLASSIFICATION_TYPE))
      stopifnot(FALSE)
    }


    save(jklm.preds, confM, bacc, file = paste0(workDir, "/", i, "_junkle_final_model_stats_preds.RData"))

    return(bacc)
  } # end of task

  return(res)
} # end of worker.f
