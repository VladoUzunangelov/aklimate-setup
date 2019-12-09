# vfriedl, July 2019 create stacked barplots for feature types used in reduced
# feature models

# run this script in the parent directory of './models'.

library(ggplot2)

this_dir_name = basename(getwd())

plot_dir <- "./"

outputfile = paste0(plot_dir, "colorfeature_proportion_barplot.png")

# colors = c('coral2', 'darkolivegreen', 'darkorchid3', 'darkgoldenrod3')


model_dir = "models/"

message("detect cutoffs for reduced models")

# get list of feature steps steps <- c(5,10,20,50,100,200,500,1000,1500)
file_list <- list.files(model_dir, pattern = "cutoff", full.names = FALSE)
steps <- sort(unique(sapply(file_list, function(x) as.numeric(strsplit(x, "_", fixed = T)[[1]][3]))))
message(steps)



message("detect feature types from full models")

# get the input file names, matching *importance.tab
file_list <- list.files(model_dir, pattern = "importance.tab", full.names = TRUE)
message(file_list)

test_dat <- read.table(paste0(model_dir, "R1:F1_aklimate_multiclass_feature_importance.tab"),
  header = TRUE, sep = "\t")
feature_types <- unique(sapply(as.character(test_dat$features), function(x) strsplit(x,
  ":", fixed = T)[[1]][2]))
print(feature_types)



message("compute datatype proportions for each reduced model cutoff")

# data frame holding the proportions for each feature type for each step
props_df <- matrix(0, nrow = length(feature_types), ncol = length(steps), dimnames = list(feature_types,
  steps))

for (step in steps) {

  # data frame to hold the counts for each feature type in each file (i.e. cv fold)
  feature_counts_df <- matrix(0, nrow = length(feature_types), ncol = length(file_list),
    dimnames = list(feature_types, file_list))

  for (file in file_list) {
    dat <- read.table(file, header = TRUE, sep = "\t")

    # if there is less features in the file than the step indicates, just take all
    # the feature in the file
    features_complete <- sapply(as.character(dat$features), function(x) strsplit(x,
      ":", fixed = T)[[1]][2])
    if (length(features_complete) >= step) {
      features <- features_complete[1:step]
    } else {
      features <- features_complete
    }

    # add counts for each feature type to the data frame
    for (ft in feature_types) {
      feature_counts_df[ft, file] <- sum(features == ft)
    }

  }

  # add proportions for each feature type to data frame
  total_num <- sum(feature_counts_df)
  for (ft in feature_types) {
    props_df[ft, toString(step)] <- sum(feature_counts_df[ft, ])/total_num
  }
}

# ensure the proportions sum up to 1
for (i in 1:length(steps)) {
  sum <- colSums(props_df)[i]
  if (sum != 1) {
    print(paste0("Warning! Proportions do not sum up to 1: ", sum))
  }
}


message("generate plot")

# sort props_df to always have the same ordering of feature types - HARDCODED FOR
# NOW
ordering <- c("GEXP", "METH", "CNVR", "MUTA")
props_df <- props_df[ordering, ]

png(outputfile)

barplot(props_df, xlab = "number of features", ylab = "proportions in top features",
  main = paste0("feature types in models\n for ", this_dir_name), legend.text = rownames(props_df),
  args.legend = list(x = "bottomright"))

# +scale_fill_manual(values = c("green", "violet", "pink"))

dev.off()
