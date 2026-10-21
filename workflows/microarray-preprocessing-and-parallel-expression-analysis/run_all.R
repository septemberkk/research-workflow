#!/usr/bin/env Rscript

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script) != 1L) stop("Use Rscript for this entry point.", call. = FALSE)
root <- dirname(normalizePath(sub("^--file=", "", script), winslash = "/", mustWork = TRUE))
for (name in c("io.R", "prepare.R", "workflow.R", "cli.R")) source(file.path(root, "src", name), local = TRUE)
tryCatch(
  array_cli(commandArgs(trailingOnly = TRUE), root),
  error = function(error) {
    cat("Error: ", conditionMessage(error), "\n", sep = "", file = stderr())
    quit(save = "no", status = 1L, runLast = FALSE)
  }
)
