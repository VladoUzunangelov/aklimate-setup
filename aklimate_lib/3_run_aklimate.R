#!/usr/bin/env Rscript
# vuzunangelov & chrisw 20190606

message("call worker on tasks")

# EDIT THIS BEFORE RUNNING !
run_tasks <- tasks[1:25]
stopifnot(length(run_tasks) > 0)
message(paste(length(run_tasks), " tasks to run"))
message(run_tasks)


# stopifnot(FALSE)

results <- worker.f(run_tasks)
names(results) <- run_tasks

results <- unlist(results)

names(results) <- as.character(run_tasks)



message("write results")

write.df(data.frame(bacc = results), "split", paste0(workDir, "/splits_balanced_accuracy.tab"))
