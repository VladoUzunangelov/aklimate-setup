#!/usr/bin/Rscript

# usage, options and doc goes here
argspec <- c("bacc_sickle_plot.R
Usage:
  Rscript bacc_sickle_plot.r balanced_accuracy_file output_file
Example:
\tRscript bacc_sickle_plot.r balanced_accuracy_file
Options:
\tbalanced_accuracy_file Columns are cutoffs for number of features allowed. Rows are repeats.
\toutput_file Name for output file. The output is a png file
")


main <- function(argv) {
  message("argv below:")
  message(argv)
  message("do something")

  if (length(argv) != 2) {
    message("ERROR: exactly two arguments required.")
    stopifnot(FALSE)
  }

  this_dir_name = basename(getwd())

  input_file = argv[1]
  outputfile = argv[2]

  message(paste0("read data for bacc_table from ", input_file))
  bacc_table = read.table(input_file, header = TRUE, row.names = 1, sep = "\t",
    check.names = FALSE)

  message("setup data")
  means = c()
  stdevs = c()
  cutoffs = colnames(bacc_table)
  x = 1:length(cutoffs)

  for (i in x) {
    means = c(means, mean(bacc_table[, i]))
    stdevs = c(stdevs, sd(bacc_table[, i]))
  }

  reps = length(rownames(bacc_table))

  message("generate plot")

  yrange = c(0, 1)

  png(outputfile)
  plot(x, means, ylim = yrange, ylab = paste0("Mean Balanced Accuracy (N=", reps,
    ")"), xaxt = "n", xlab = "Number of Features", main = paste0("# Features vs. Balanced Accuracy\nfor ",
    this_dir_name))
  grid()
  lines(x, means, col = "red")
  axis(1, at = x, labels = colnames(bacc_table))
  arrows(x, means - stdevs, x, means + stdevs, length = 0.05, angle = 90, code = 3)
  dev.off()

}

main(commandArgs(TRUE))
