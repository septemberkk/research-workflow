mirna_help <- function() {
  paste(c(
    "miRNA Differential Expression & Target Enrichment",
    "Independent two-group limma analysis; reviewed target-table handoff, no enrichment calculation.",
    "",
    "Usage:",
    "  Rscript run_all.R --help",
    "  Rscript run_all.R --check-environment",
    "  Rscript run_all.R --example --output OUTPUT [--plots]",
    "  Rscript run_all.R --expression FILE --samples FILE --reference GROUP --comparison GROUP --input-scale log2 --output OUTPUT [options]",
    "",
    "Options:",
    "  --normalization METHOD  quantile (default) or none for already normalized log2 values.",
    "  --duplicates POLICY     error (default) or mean for repeated miRNA IDs.",
    "  --fdr NUMBER            BH value strictly below this threshold (0.05).",
    "  --logfc NUMBER          Absolute log2 fold change strictly above this threshold (1).",
    "  --targets FILE          Optional reviewed target-evidence TSV.",
    "  --plots                 Base R volcano and row-z-score heatmap.",
    "  --heatmap-rows N        Maximum selected heatmap rows (40).",
    "",
    "Expression: miRNA_id followed by numeric sample columns.",
    "Samples: sample_id and group, with exactly two independent groups.",
    "Output must be new; its parent must exist. No existing files are overwritten.",
    "No raw counts, background correction, paired/batch model, database queries, package installation or network access.",
    "Input locations are supplied at runtime and are not recorded in output metadata."
  ), collapse = "\n")
}

mirna_parse <- function(args) {
  flags <- c("--help", "--check-environment", "--example", "--plots")
  values <- c("--expression", "--samples", "--reference", "--comparison", "--input-scale",
               "--output", "--normalization", "--duplicates", "--fdr", "--logfc",
               "--targets", "--heatmap-rows")
  parsed <- list()
  i <- 1L
  while (i <= length(args)) {
    option <- args[i]
    if (!option %in% c(flags, values)) mirna_fail("Unknown option. See --help.")
    key <- sub("^--", "", option)
    if (key %in% names(parsed)) mirna_fail("Duplicate option.")
    if (option %in% flags) parsed[[key]] <- TRUE else {
      if (i == length(args) || grepl("^--", args[i + 1L])) mirna_fail("Missing option value.")
      i <- i + 1L
      parsed[[key]] <- mirna_text(args[i], key)
    }
    i <- i + 1L
  }
  for (action in c("help", "check-environment")) {
    if (isTRUE(parsed[[action]])) {
      if (length(parsed) != 1L) mirna_fail("Help and environment actions must be used alone.")
      return(list(action = action))
    }
  }
  if (isTRUE(parsed$example)) {
    if (any(c("expression", "samples", "reference", "comparison", "input-scale", "targets") %in% names(parsed))) {
      mirna_fail("Example mode cannot be combined with user inputs or group labels.")
    }
    parsed$reference <- "Control"
    parsed$comparison <- "Treatment"
    parsed[["input-scale"]] <- "log2"
  }
  required <- c(if (!isTRUE(parsed$example)) c("expression", "samples"),
                 "reference", "comparison", "input-scale", "output")
  if (any(!required %in% names(parsed))) mirna_fail("Required options are missing. See --help.")
  number <- function(key, default) {
    if (is.null(parsed[[key]])) default else suppressWarnings(as.numeric(parsed[[key]]))
  }
  parsed$config <- mirna_config(
    parsed$reference, parsed$comparison, parsed[["input-scale"]],
    normalization = if (is.null(parsed$normalization)) "quantile" else parsed$normalization,
    duplicates = if (is.null(parsed$duplicates)) "error" else parsed$duplicates,
    fdr = number("fdr", .05), logfc = number("logfc", 1),
    heatmap_rows = number("heatmap-rows", 40)
  )
  parsed$action <- "run"
  parsed$example <- isTRUE(parsed$example)
  parsed$plots <- isTRUE(parsed$plots)
  parsed
}

mirna_cli <- function(args, root) {
  options <- mirna_parse(args)
  if (options$action == "help") return(cat(mirna_help(), "\n", sep = ""))
  if (options$action == "check-environment") {
    print(mirna_environment(), row.names = FALSE)
    mirna_require_limma()
    return(invisible(NULL))
  }
  expression <- if (options$example) file.path(root, "examples", "expression.tsv") else options$expression
  samples <- if (options$example) file.path(root, "examples", "samples.tsv") else options$samples
  targets <- if (options$example) file.path(root, "examples", "reviewed_targets.tsv") else options$targets
  result <- mirna_run(expression, samples, options$output, options$config, root,
                       targets_file = targets, plots = options$plots)
  cat("Differential expression complete: ", nrow(result$statistics), " features; ",
       sum(result$statistics$selected), " selected. Target query and functional enrichment: NOT RUN.\n", sep = "")
  invisible(result)
}
