#!/usr/bin/Rscript

# Get the sample label predictions from AKLIMATE output file.

# usage, options and doc goes here
argspec <- c("get_sample_predictions.R
Usage:
    Rscript get_sample_predictions.R filename
Example:
\tRscript get_sample_predictions.R ABCD_junkle_final_model_stats_preds.RData
Options:
\tfilename is an AKLIMATE output file with name like '_junkle_final_model_stats_preds.RData'. It contains an object jklm.pred.
")

main <- function(argv) {
    load(argv[1])
    write.df(data.frame(as.list(jklm.preds)),row.names.id = "sampleID", "preds.tmp")
}

main(commandArgs(TRUE))
