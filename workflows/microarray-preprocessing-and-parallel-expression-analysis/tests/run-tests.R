#!/usr/bin/env Rscript

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script) != 1L) stop("Use Rscript for this entry point.", call. = FALSE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", script), winslash = "/", mustWork = TRUE)))
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) args <- "--unit"
if (length(args) != 1L || !args %in% c("--unit", "--integration", "--all")) {
  stop("Use --unit, --integration or --all.", call. = FALSE)
}
for (name in c("io.R", "prepare.R", "workflow.R", "cli.R")) source(file.path(root, "src", name), local = TRUE)
source(file.path(root, "tests", "helpers.R"), local = TRUE)
passed <- 0L
failed <- 0L
test <- function(name, code) {
  tryCatch({
    force(code)
    passed <<- passed + 1L
    cat("PASS ", name, "\n", sep = "")
  }, error = function(error) {
    failed <<- failed + 1L
    cat("FAIL ", name, ": ", conditionMessage(error), "\n", sep = "")
  })
}
expect <- function(value) {
  if (!isTRUE(value)) stop("Expected TRUE.", call. = FALSE)
}
expect_equal <- function(actual, expected, tolerance = 1e-10) {
  result <- all.equal(actual, expected, tolerance = tolerance, check.attributes = TRUE)
  if (!isTRUE(result)) stop(paste(result, collapse = "; "), call. = FALSE)
}
expect_error <- function(code, pattern) {
  error <- tryCatch({force(code); NULL}, error = identity)
  if (is.null(error) || !grepl(pattern, conditionMessage(error), fixed = TRUE)) {
    stop("Expected error with the specified message.", call. = FALSE)
  }
}
if (args %in% c("--unit", "--all")) source(file.path(root, "tests", "unit.R"), local = TRUE)
blocked <- FALSE
if (args %in% c("--integration", "--all")) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    blocked <- TRUE
    cat("BLOCKED: quantile-engine checks were not executed. No packages installed.\n")
  } else source(file.path(root, "tests", "integration.R"), local = TRUE)
}
cat("SUMMARY passed=", passed, " failed=", failed, " integration_blocked=", blocked, "\n", sep = "")
quit(save = "no", status = if (failed) 1L else if (blocked) 2L else 0L, runLast = FALSE)
