# View(confM$byClass) View(confM$table) (Cohort, Attributes, Readily, Identified,
# Now, Algorithm) Gabe & Jackie & Verrena & Chris
# ****************************************************************************************************************
# This script continues the process of the 'find_balance_accuracy_cohorts.R'
# script. Instead of looking at the balance accuracies of each cohort however, we
# find the balance accuracy of the invidivual subcorts in each cohort and
# individually plot each tumor types sub tumor types.
# ****************************************************************************************************************
# file_name_list paths: in ucsc
# terminal:'/scratch/for_gchavez/aklimate_results/',tt,'/models/',sep='' in gabes
# computer: /Users/user/Desktop/BD2K_project/data/ in jackies computer:
# /Users/jacquelynroger/Documents/research/RMI/gabe/data/ tumor_type_list: the
# types of cohorts we're comparing stats = confM[[4]] bal_accs = stats[,11] brca
# classes path to brca

reduced_model_cutoff <- 50

# outputfile <- paste0('/scratch/for_gchavez/aklimate_results', 'sub_', tt,
# '-plot.pdf')
outputfile <- paste0("bal_acc_subtype_50_cutoff.png")

subtype_name_mapping_file <- "/scratch/for_gchavez/aklimate_results/subtypes_mapping.tsv"

this_dir_name = basename(getwd())


tumor_type_list = c("brca", "coadread", "lgggbm", "thym", "ucec")
colors = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D")
df_cohort = c()
class_types = c("bal_accs_class1", "bal_accs_class2", "bal_accs_class3", "bal_accs_class4")



tt = this_dir_name
model_dir = "models/"

# filepath = paste('/scratch/for_gchavez/aklimate_results/', tt, '/models/', sep
# = '') # file_name_list: Takes all the files and lists them file_name_list =
# list.files(path = filepath)

file_name_list <- list.files(model_dir, pattern = paste0("_cutoff_", reduced_model_cutoff,
  "_rf_reduced_model_stats.RData"), full.names = FALSE)

# num_files: The number of files we have
num_files = length(file_name_list)
message(paste0("num_files: ", num_files))
message(file_name_list)


# Initialize empty list to store the aggregate of all of the balance accuracies
# in the cohort subtypes
crossval_bal_acc = c()
for (i in 1:num_files) {
  file_name = file_name_list[(i)]
  # LOOP B: This if statement makes sure we only go into the files that are
  # specifically AKLIMATE output files

  load(paste0(model_dir, file_name))
  message(paste0("loaded file ", file_name))

  # Initialize empty list to store the aggregate of all of the balance accuracies
  # in this cohort stats = confM[[4]]
  if (length(rownames(confM[["byClass"]])) == 1) {
    message("detected results for binary classification model")
    stats <- unname(confM[["byClass"]]["Balanced Accuracy"])

  } else {
    message("detected results for multiclass model")
    stats <- (confM$byClass)[, 11]
  }
  class_names <- names(stats)
  message(stats)
  message(class_names)
  bal_accs_class = c()
  for (w in 1:length(class_names)) {
    class_name <- class_names[w]
    message(class_name)
    bal_accs = stats[w]
    message(bal_accs)
    bal_accs_class = c(bal_accs_class, bal_accs)
  }
  crossval_bal_acc = rbind(crossval_bal_acc, bal_accs_class)

}

message(crossval_bal_acc)
message(colnames(crossval_bal_acc))

message("got bal acc stats")

mean_data = c()
std_data = c()
z = 1:ncol(crossval_bal_acc)
for (i in z) {
  mean_data = c(mean_data, mean(crossval_bal_acc[, i]))
  std_data = c(std_data, sd(crossval_bal_acc[, i]))
}

message(mean_data)
message(std_data)

if (FALSE) {

  message("map class name to subtype name")
  # TODO this part needs fixing this lists the first three letters in each
  # subtype_number substr(dat$subtype_number, 1, 3)
  sub_cohort_names = c()
  dat = read.table(subtype_name_mapping_file, header = TRUE, sep = "\t", check.names = FALSE,
    stringsAsFactors = FALSE)
  for (i in 1:length(dat$subtype_number)) {
    if (toupper(substr(tt, 1, 3)) == substr(dat$subtype_number[i], 1, 3)) {
      # print(dat[i,])
      sub_cohort_names = c(sub_cohort_names, dat$subtype_name[i])
    }
  }

}

png(outputfile)
text_main <- paste0("Balanced accuracy by subtype\nfor ", tt, "\n(", reduced_model_cutoff,
  " features)")
# yrange <- c(-0.35,1.1)
yrange <- c(0, 1.06)
x = barplot(colMeans(crossval_bal_acc), main = text_main, xlab = "Subtype", ylab = "Balanced Accuracy",
  ylim = yrange, col = colors, border = "white", yaxt = "n")

# this adds error bars
arrows(x, mean_data - std_data, x, mean_data + std_data, length = 0.05, angle = 90,
  code = 3)

# this round the BalAcc for each bar to the 2nd decimal digit
y = round(colMeans(crossval_bal_acc), digits = 2)

# this creates the y axis
axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = as.character(c(0, 0.2, 0.4, 0.6,
  0.8, 1)))

# this moves the text for the BalAcc on the bars
text(x - 0.3, colMeans(crossval_bal_acc) + 0.03, labels = as.character(y))


# TODO this part needs fixing this is the text for the subcohorts
if (tt == "lgggbm" | tt == "ucec") {
  text(x - 0.25, rep(-0.15, length.out = length(x)), labels = sub_cohort_names[1:4],
    srt = 60)
}
if (tt == "brca" | tt == "thym") {
  text(x - 0.25, rep(-0.1, length.out = length(x)), labels = sub_cohort_names[1:4],
    srt = 60)
}
if (tt == "coadread") {
  text(x - 0.25, rep(-0.1, length.out = length(x)), labels = sub_cohort_names[1:4],
    srt = 60, cex = 0.8)
}

dev.off()



# how to run cd /scratch/for_gchavez/aklimate_results/lib/AKLIMATE_scripts
# Rscript sub_bal_acc_v1.R
