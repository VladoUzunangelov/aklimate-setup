#!/usr/bin/env Rscript
# vuzunangelov & chrisw 20190606

message("call worker on tasks")

# EDIT THIS BEFORE RUNNING !
run_tasks <- tasks[1:25]
stopifnot(length(run_tasks) > 0)
message(paste(length(run_tasks), " tasks to run"))
message(run_tasks)


# stopifnot(FALSE)


message("==> train/test models")
t0 <- Sys.time()

results <- worker.f(run_tasks)

t1 <- Sys.time()
message("==> got results object")
dt <- difftime(t1, t0, units = "secs")
message(paste0("experiment run time (secs): ", dt))


names(results) <- run_tasks

results <- unlist(results)

names(results) <- as.character(run_tasks)



message("write results")

write.df(data.frame(bacc = results), "split", paste0(workDir, "/splits_balanced_accuracy.tab"))
