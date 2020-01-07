#!/usr/bin/Rscript

main <- function(argv) {
    load(argv[1])
	cat(sprintf("%f\n", bacc))
}

main(commandArgs(TRUE))
