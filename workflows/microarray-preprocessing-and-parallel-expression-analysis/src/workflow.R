array_environment <- function() {
  available <- requireNamespace("limma", quietly = TRUE)
  data.frame(package = "limma", available = available,
              role = "quantile normalization only",
              version = if (available) as.character(utils::packageVersion("limma")) else "unavailable")
}

array_git <- function(root) {
  git <- Sys.which("git")
  if (!nzchar(git)) return("unavailable")
  result <- suppressWarnings(system2(git, c("-C", shQuote(root), "rev-parse", "HEAD"),
                                      stdout = TRUE, stderr = FALSE))
  if (length(result) == 1L && grepl("^[0-9a-f]{40}$", result)) result else "unavailable"
}

array_run <- function(files, output, config, root) {
  started <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  array_output(output)
  roles <- c("mrna_expression", "mrna_annotation", "mrna_samples",
              "mirna_expression", "mirna_annotation", "mirna_samples")
  if (!is.character(files) || !setequal(names(files), roles) || anyDuplicated(names(files))) {
    array_fail("Exactly six named expression, annotation and sample inputs are required.")
  }
  files <- files[roles]
  for (file in files) {
    if (!file.exists(file) || dir.exists(file)) array_fail("A required input is missing.")
  }
  manifests <- array_manifest(files, roles)
  results <- list()
  for (layer in c("mrna", "mirna")) {
    result <- array_prepare(
      array_read(files[[paste0(layer, "_expression")]], "expression"),
      array_read(files[[paste0(layer, "_annotation")]], "annotation"),
      array_read(files[[paste0(layer, "_samples")]], "samples"), config[[layer]]$scale
    )
    result$processed <- array_normalize(result$collapsed, config[[layer]]$normalization)
    result$distribution <- rbind(array_distribution(result$collapsed, "collapsed"),
                                  array_distribution(result$processed, "processed"))
    results[[layer]] <- result
  }
  if (!identical(manifests, array_manifest(files, roles))) array_fail("Inputs changed during processing; no output was created.")
  source_roles <- c("run_all.R", paste0("src/", c("io.R", "prepare.R", "workflow.R", "cli.R")))
  source_manifest <- array_manifest(file.path(root, source_roles), source_roles)
  if (!dir.create(output, recursive = FALSE)) array_fail("Cannot create a new output directory.")
  status <- "failed"
  on.exit({
    metadata <- c(workflow = "microarray-preprocessing-and-parallel-expression-analysis",
                   status = status, started_utc = started,
                   finished_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
                   R_version = as.character(getRversion()), git_commit = array_git(root),
                   limma_version = array_environment()$version,
                   mrna_scale = config$mrna$scale, mirna_scale = config$mirna$scale,
                   mrna_normalization = config$mrna$normalization,
                   mirna_normalization = config$mirna$normalization,
                   collapsed_on = "declared_supplied_scale", log_transform = "not_run",
                   background_correction = "not_run", differential_expression = "not_run",
                   cross_layer_pairing = "not_inferred", regulatory_interactions = "not_inferred",
                   cross_layer_integration = "not_run", plots_requested = config$plots)
    array_write(data.frame(key = names(metadata), value = unname(metadata)), file.path(output, "run_metadata.tsv"))
    exports <- list.files(output, full.names = TRUE)
    exports <- exports[basename(exports) != "output_manifest.tsv"]
    array_write(array_manifest(exports, basename(exports)), file.path(output, "output_manifest.tsv"))
  })
  array_write(manifests, file.path(output, "input_manifest.tsv"))
  array_write(source_manifest, file.path(output, "source_manifest.tsv"))
  for (layer in names(results)) {
    result <- results[[layer]]
    write <- function(table, suffix) array_write(table, file.path(output, paste0(layer, "_", suffix, ".tsv")))
    write(result$mapping, "probe_mapping")
    write(result$members, "feature_members")
    write(result$counts, "feature_counts")
    write(result$samples, "samples_supplied")
    write(result$selected, "samples_selected")
    write(result$annotation, "annotation_supplied")
    write(result$summary, "preparation_summary")
    write(result$distribution, "sample_distributions")
    for (stage in c("collapsed", "processed")) {
      matrix <- result[[stage]]
      write(data.frame(feature_id = rownames(matrix), matrix, row.names = NULL, check.names = FALSE),
             paste0(stage, "_expression"))
    }
    if (config$plots) array_plot(layer, result, config[[layer]]$scale,
                                  file.path(output, paste0(layer, "_distributions.pdf")))
  }
  array_write(data.frame(layer = names(results), status = if (config$plots) "completed" else "not_requested"),
                file.path(output, "plot_status.tsv"))
  status <- "completed"
  invisible(results)
}
