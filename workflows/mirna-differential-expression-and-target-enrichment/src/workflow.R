mirna_environment <- function() {
  data.frame(package = "limma", available = requireNamespace("limma", quietly = TRUE),
              version = if (requireNamespace("limma", quietly = TRUE)) {
                as.character(utils::packageVersion("limma"))
              } else "unavailable")
}

mirna_require_limma <- function() {
  if (!requireNamespace("limma", quietly = TRUE)) {
    mirna_fail("limma is unavailable. No packages were installed and no analysis was run.")
  }
}

mirna_fit_limma <- function(inputs, config) {
  mirna_require_limma()
  normalized <- if (config$normalization == "quantile") {
    limma::normalizeBetweenArrays(inputs$expression, method = "quantile")
  } else inputs$expression
  fit <- limma::lmFit(normalized, inputs$design)
  contrast <- matrix(c(-1, 1), ncol = 1L,
                      dimnames = list(colnames(inputs$design), "comparison_minus_reference"))
  fit <- limma::contrasts.fit(fit, contrast)
  fit <- limma::eBayes(fit, trend = FALSE, robust = FALSE)
  table <- limma::topTable(fit, coef = 1, number = Inf, adjust.method = "BH", sort.by = "P")
  list(expression = normalized, statistics = table)
}

mirna_statistics <- function(table, ids, config) {
  required <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
  if (!is.data.frame(table) || !all(required %in% names(table)) ||
      nrow(table) != length(ids) || anyDuplicated(rownames(table)) ||
      !setequal(rownames(table), ids)) {
    mirna_fail("limma statistics do not cover the complete assayed feature set.")
  }
  if (any(!vapply(table[required], is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(table[required]))) ||
      any(table$P.Value < 0 | table$P.Value > 1 | table$adj.P.Val < 0 | table$adj.P.Val > 1)) {
    mirna_fail("limma returned invalid or nonfinite statistics.")
  }
  result <- data.frame(miRNA_id = rownames(table), table[required],
                        row.names = NULL, check.names = FALSE)
  result <- result[order(result$P.Value, result$miRNA_id, method = "radix"), , drop = FALSE]
  rownames(result) <- NULL
  result$selected <- result$adj.P.Val < config$fdr & abs(result$logFC) > config$logfc
  result$direction <- ifelse(result$selected, ifelse(result$logFC > 0, "up", "down"), "not_selected")
  result
}

mirna_manifest <- function(files, roles) {
  data.frame(role = roles, bytes = unname(file.info(files)$size),
              md5 = unname(tools::md5sum(files)), stringsAsFactors = FALSE)
}

mirna_source_manifest <- function(root) {
  files <- c("run_all.R", file.path("src", c("input.R", "targets.R", "plots.R", "workflow.R", "cli.R")))
  result <- mirna_manifest(file.path(root, files), gsub("\\\\", "/", files))
  result
}

mirna_git_commit <- function(root) {
  git <- Sys.which("git")
  if (!nzchar(git)) return("unavailable")
  result <- suppressWarnings(system2(git, c("-C", shQuote(root), "rev-parse", "HEAD"),
                                      stdout = TRUE, stderr = FALSE))
  if (length(result) != 1L || !grepl("^[0-9a-f]{40}$", result)) return("unavailable")
  result
}

mirna_run <- function(expression_file, samples_file, output, config, root,
                      targets_file = NULL, plots = FALSE) {
  started <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  mirna_new_output(output)
  input_files <- c(expression_file, samples_file, if (!is.null(targets_file)) targets_file)
  input_roles <- c("expression", "samples", if (!is.null(targets_file)) "reviewed_targets")
  for (file in input_files) {
    if (!file.exists(file) || dir.exists(file)) mirna_fail("A required input file is missing.")
  }
  input_manifest <- mirna_manifest(input_files, input_roles)
  expression <- mirna_read_tsv(expression_file, "expression")
  samples <- mirna_read_tsv(samples_file, "samples")
  inputs <- mirna_inputs(expression, samples, config)
  targets <- if (is.null(targets_file)) NULL else {
    mirna_targets(mirna_read_tsv(targets_file, "target evidence"), rownames(inputs$expression))
  }
  mirna_require_limma()
  fitted <- mirna_fit_limma(inputs, config)
  statistics <- mirna_statistics(fitted$statistics, rownames(inputs$expression), config)
  target_exports <- mirna_target_exports(targets, statistics)
  source_manifest <- mirna_source_manifest(root)
  if (!identical(input_manifest, mirna_manifest(input_files, input_roles))) {
    mirna_fail("Input files changed during analysis; no output was created.")
  }
  if (!dir.create(output, recursive = FALSE)) mirna_fail("Cannot create a new output directory.")
  status <- "failed"
  on.exit({
    metadata <- c(workflow = "mirna-differential-expression-and-target-enrichment",
                   status = status, finished_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
                   started_utc = started, R_version = as.character(getRversion()),
                   limma_version = as.character(utils::packageVersion("limma")),
                   git_commit = mirna_git_commit(root), input_scale = config$input_scale,
                   normalization = config$normalization, duplicate_policy = config$duplicates,
                   reference = config$reference, comparison = config$comparison,
                   reference_samples = sum(samples$group == config$reference),
                   comparison_samples = sum(samples$group == config$comparison),
                   input_rows = inputs$input_rows, collapsed_rows = inputs$collapsed_rows,
                   tested_features = nrow(statistics), selected_features = sum(statistics$selected),
                   fdr_strictly_below = config$fdr, absolute_logfc_strictly_above = config$logfc,
                   plots_requested = plots, heatmap_rows = config$heatmap_rows,
                   target_evidence = if (is.null(targets)) "not_supplied" else "reviewed_table",
                   target_query = "not_run", functional_enrichment = "not_run",
                   enrichment_background = "not_defined", target_expression_direction = "not_inferred")
    mirna_write_tsv(data.frame(key = names(metadata), value = unname(metadata)),
                     file.path(output, "run_metadata.tsv"))
    files <- list.files(output, full.names = TRUE, recursive = FALSE)
    files <- files[basename(files) != "output_manifest.tsv"]
    mirna_write_tsv(mirna_manifest(files, basename(files)), file.path(output, "output_manifest.tsv"))
  })
  mirna_write_tsv(input_manifest, file.path(output, "input_manifest.tsv"))
  mirna_write_tsv(source_manifest, file.path(output, "source_manifest.tsv"))
  mirna_write_tsv(samples, file.path(output, "samples_used.tsv"))
  mirna_write_tsv(data.frame(sample_id = rownames(inputs$design), inputs$design,
                             row.names = NULL, check.names = FALSE),
                   file.path(output, "design.tsv"))
  mirna_write_tsv(data.frame(miRNA_id = rownames(fitted$expression), fitted$expression,
                             row.names = NULL, check.names = FALSE),
                   file.path(output, "normalized_expression.tsv"))
  mirna_write_tsv(statistics, file.path(output, "differential_all.tsv"))
  mirna_write_tsv(statistics[statistics$selected, , drop = FALSE],
                   file.path(output, "differential_selected.tsv"))
  mirna_write_ids(statistics$miRNA_id[statistics$direction == "up"], file.path(output, "miRNAs_up.txt"))
  mirna_write_ids(statistics$miRNA_id[statistics$direction == "down"], file.path(output, "miRNAs_down.txt"))
  if (!is.null(target_exports)) {
    mirna_write_tsv(target_exports$review, file.path(output, "target_review.tsv"))
    mirna_write_tsv(target_exports$evidence, file.path(output, "target_evidence_selected.tsv"))
    mirna_write_ids(target_exports$all, file.path(output, "target_genes_all.txt"))
    mirna_write_ids(target_exports$up, file.path(output, "targets_of_up_miRNAs.txt"))
    mirna_write_ids(target_exports$down, file.path(output, "targets_of_down_miRNAs.txt"))
  }
  plot_status <- mirna_plots(fitted$expression, statistics, config, output, plots)
  mirna_write_tsv(plot_status, file.path(output, "plot_status.tsv"))
  status <- "completed"
  invisible(list(statistics = statistics, expression = fitted$expression,
                  targets = target_exports, plots = plot_status))
}
