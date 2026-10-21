#!/usr/bin/env Rscript

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script) != 1L) stop("Use Rscript to run this entry point.", call. = FALSE)
root <- dirname(normalizePath(sub("^--file=", "", script), winslash = "/", mustWork = TRUE))
for (file in c("input.R", "targets.R", "plots.R", "workflow.R", "cli.R")) {
  source(file.path(root, "src", file), local = TRUE)
}
tryCatch(
  mirna_cli(commandArgs(trailingOnly = TRUE), root),
  error = function(error) {
    cat("Error: ", conditionMessage(error), "\n", sep = "", file = stderr())
    quit(save = "no", status = 1L, runLast = FALSE)
  }
)
