test("workflow equals a direct limma reference analysis with quantile normalization", {
  mirna_temp(function(directory) {
    fixture <- fixture_inputs()
    inputs <- mirna_inputs(fixture$expression, fixture$samples, fixture_config())
    normalized <- limma::normalizeBetweenArrays(inputs$expression, method = "quantile")
    fit <- limma::lmFit(normalized, inputs$design)
    contrast <- matrix(c(-1, 1), ncol = 1L,
                        dimnames = list(colnames(inputs$design), "comparison_minus_reference"))
    fit <- limma::eBayes(limma::contrasts.fit(fit, contrast), trend = FALSE, robust = FALSE)
    reference <- limma::topTable(fit, coef = 1, number = Inf, adjust.method = "BH", sort.by = "P")
    expected <- mirna_statistics(reference, rownames(normalized), fixture_config())
    output <- file.path(directory, "result")
    result <- mirna_run(file.path(root, "examples", "expression.tsv"),
                         file.path(root, "examples", "samples.tsv"), output, fixture_config(), root,
                         file.path(root, "examples", "reviewed_targets.tsv"), plots = TRUE)
    expect_equal(result$statistics, expected)
    expect_equal(result$expression, normalized)
    expect(sum(result$statistics$selected) > 0)
    expect(all(result$plots$status == "completed"))
    metadata <- read.delim(file.path(output, "run_metadata.tsv"), check.names = FALSE)
    expect_equal(metadata$value[metadata$key == "status"], "completed")
    expect_equal(metadata$value[metadata$key == "functional_enrichment"], "not_run")
    expect_equal(metadata$value[metadata$key == "enrichment_background"], "not_defined")
    expect_equal(metadata$value[metadata$key == "target_expression_direction"], "not_inferred")
    manifest <- read.delim(file.path(output, "output_manifest.tsv"), check.names = FALSE)
    expect_equal(unname(tools::md5sum(file.path(output, manifest$role))), manifest$md5)
    for (name in c("input_manifest.tsv", "run_metadata.tsv")) {
      expect(!any(grepl(root, readLines(file.path(output, name)), fixed = TRUE)))
      expect(!any(grepl(directory, readLines(file.path(output, name)), fixed = TRUE)))
    }
  })
})
test("no-normalization mode preserves input values and exact contrast direction", {
  fixture <- fixture_inputs()
  config <- fixture_config(normalization = "none")
  inputs <- mirna_inputs(fixture$expression, fixture$samples, config)
  result <- mirna_fit_limma(inputs, config)
  expect_equal(result$expression, inputs$expression)
  expected <- rowMeans(inputs$expression[, inputs$samples$group == "Treatment", drop = FALSE]) -
    rowMeans(inputs$expression[, inputs$samples$group == "Control", drop = FALSE])
  expect_equal(unname(result$statistics[match(names(expected), rownames(result$statistics)), "logFC"]),
                unname(expected))
})
test("empty selections export valid headers and no invented targets or heatmap", {
  mirna_temp(function(directory) {
    output <- file.path(directory, "empty")
    result <- mirna_run(file.path(root, "examples", "expression.tsv"),
                         file.path(root, "examples", "samples.tsv"), output,
                         fixture_config(logfc = 100), root, plots = TRUE)
    expect_equal(sum(result$statistics$selected), 0L)
    expect_equal(length(readLines(file.path(output, "differential_selected.tsv"))), 1L)
    expect_equal(file.info(file.path(output, "miRNAs_up.txt"))$size, 0)
    expect_equal(file.info(file.path(output, "miRNAs_down.txt"))$size, 0)
    expect(!file.exists(file.path(output, "target_evidence_selected.tsv")))
    expect(!file.exists(file.path(output, "heatmap.pdf")))
    expect_equal(result$plots$status[2], "skipped_no_selected_variable_rows")
  })
})
test("input byte manifests remain stable and overwriting remains impossible", {
  mirna_temp(function(directory) {
    files <- file.path(root, "examples", c("expression.tsv", "samples.tsv", "reviewed_targets.tsv"))
    before <- tools::md5sum(files)
    output <- file.path(directory, "result")
    mirna_run(files[1], files[2], output, fixture_config(), root, files[3])
    expect_equal(tools::md5sum(files), before)
    metadata <- readLines(file.path(output, "run_metadata.tsv"))
    expect_error(mirna_run(files[1], files[2], output, fixture_config(), root), "already exists")
    expect_equal(readLines(file.path(output, "run_metadata.tsv")), metadata)
  })
})
