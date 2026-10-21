array_help <- function() {
  paste(c(
    "Microarray Preprocessing & Parallel Expression Analysis",
    "Independent mRNA and miRNA preparation; no differential-expression model or cross-layer inference.",
    "",
    "Usage:",
    "  Rscript run_all.R --help",
    "  Rscript run_all.R --check-environment",
    "  Rscript run_all.R --example --output OUTPUT [--plots]",
    "",
    "User inputs require all six --mrna-expression, --mrna-annotation, --mrna-samples,",
    "--mirna-expression, --mirna-annotation and --mirna-samples FILE options,",
    "plus --mrna-scale and --mirna-scale (intensity or log2), and --output OUTPUT.",
    "",
    "Options:",
    "  --mrna-normalization METHOD   none (default) or quantile.",
    "  --mirna-normalization METHOD  quantile (default) or none.",
    "  --plots                       Independent per-layer expression distributions.",
    "",
    "Expression headers: probe_id and sample IDs.",
    "Annotation headers: probe_id and feature_id; blank/ambiguous mappings are reported.",
    "Samples headers: sample_id, group, include (yes/no).",
    "All input locations are runtime arguments and are excluded from exported metadata.",
    "Output must be new, with an existing parent. No installation, download or input rewriting."
  ), collapse = "\n")
}

array_parse <- function(args) {
  flags <- c("--help", "--check-environment", "--example", "--plots")
  input_keys <- c("mrna-expression", "mrna-annotation", "mrna-samples",
                   "mirna-expression", "mirna-annotation", "mirna-samples")
  values <- paste0("--", c(input_keys, "mrna-scale", "mirna-scale",
                            "mrna-normalization", "mirna-normalization", "output"))
  options <- list()
  index <- 1L
  while (index <= length(args)) {
    key <- args[index]
    if (!key %in% c(flags, values)) array_fail("Unknown option. See --help.")
    name <- sub("^--", "", key)
    if (name %in% names(options)) array_fail("Duplicate option.")
    if (key %in% flags) options[[name]] <- TRUE else {
      if (index == length(args) || grepl("^--", args[index + 1L])) array_fail("Missing option value.")
      index <- index + 1L
      options[[name]] <- array_text(args[index], name)
    }
    index <- index + 1L
  }
  for (action in c("help", "check-environment")) {
    if (isTRUE(options[[action]])) {
      if (length(options) != 1L) array_fail("Help and environment actions must stand alone.")
      return(list(action = action))
    }
  }
  if (isTRUE(options$example)) {
    if (any(c(input_keys, "mrna-scale", "mirna-scale") %in% names(options))) {
      array_fail("Example mode cannot be combined with user inputs or scales.")
    }
    options[["mrna-scale"]] <- "intensity"
    options[["mirna-scale"]] <- "intensity"
  }
  required <- c(if (!isTRUE(options$example)) input_keys, "mrna-scale", "mirna-scale", "output")
  if (any(!required %in% names(options))) array_fail("Required options are missing. See --help.")
  options$config <- array_config(options[["mrna-scale"]], options[["mirna-scale"]],
                                  if (is.null(options[["mrna-normalization"]])) "none" else options[["mrna-normalization"]],
                                  if (is.null(options[["mirna-normalization"]])) "quantile" else options[["mirna-normalization"]],
                                  plots = isTRUE(options$plots))
  options$action <- "run"
  options$example <- isTRUE(options$example)
  options
}

array_example_files <- function(root) {
  roles <- c("mrna_expression", "mrna_annotation", "mrna_samples",
              "mirna_expression", "mirna_annotation", "mirna_samples")
  stats::setNames(file.path(root, "examples", paste0(roles, ".tsv")), roles)
}

array_cli <- function(args, root) {
  options <- array_parse(args)
  if (options$action == "help") return(cat(array_help(), "\n", sep = ""))
  if (options$action == "check-environment") {
    print(array_environment(), row.names = FALSE)
    return(invisible(NULL))
  }
  roles <- c("mrna_expression", "mrna_annotation", "mrna_samples",
              "mirna_expression", "mirna_annotation", "mirna_samples")
  files <- if (options$example) array_example_files(root) else {
    stats::setNames(vapply(gsub("_", "-", roles), function(key) options[[key]], character(1)), roles)
  }
  result <- array_run(files, options$output, options$config, root)
  cat("Independent mRNA and miRNA preparation complete. Differential expression and cross-layer inference: NOT RUN.\n")
  invisible(result)
}
