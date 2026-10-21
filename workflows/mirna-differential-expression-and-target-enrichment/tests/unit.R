test("default configuration states scale, direction and strict screen", {
  config <- fixture_config()
  expect_equal(config$normalization, "quantile")
  expect_equal(config$duplicates, "error")
  expect_equal(config$fdr, .05)
  expect_equal(config$logfc, 1)
})
for (value in c("raw", "counts", "auto", "ln", "LOG2", "")) {
  test(paste("reject undeclared or unsupported scale", value), {
    expect_error(mirna_config("Control", "Treatment", value), "log2")
  })
}
for (arguments in list(list(fdr = 0), list(fdr = 2), list(fdr = NA_real_),
                       list(logfc = -1), list(logfc = Inf), list(heatmap_rows = 0),
                       list(heatmap_rows = 1.5), list(heatmap_rows = 1001))) {
  test(paste("reject invalid numeric setting", names(arguments)[1]), {
    expect_error(do.call(fixture_config, arguments), "Invalid")
  })
}
test("same, missing and unsafe group labels are rejected", {
  expect_error(mirna_config("A", "A", "log2"), "different")
  expect_error(mirna_config(NULL, "A", "log2"), "different")
  expect_error(mirna_config("A ", "B", "log2"), "boundary")
  expect_error(mirna_config("=A", "B", "log2"), "formula-like")
})
test("unsupported normalization and duplicate policies are rejected", {
  expect_error(fixture_config(normalization = "auto"), "Normalization")
  expect_error(fixture_config(duplicates = "first"), "Duplicate policy")
})
for (value in c("", " A", "A ", "=formula", "+formula", "@formula", "-formula", "A\u200b", "A\nB")) {
  test("identifier and cell safety", expect_error(mirna_safe_cells(value, "identifier"),
                                                   if (grepl("^[=+@-]", value) || grepl("\u200b", value)) "formula-like" else "text"))
}
test("literal numeric-looking identifiers are not coerced", {
  mirna_safe_cells(c("001", "NA", "MIR_001"), "identifier")
})
for (value in c('"=payload"', '""', 'MIR_"A"')) {
  test("quoted fields cannot bypass spreadsheet safety", {
    expect_error(mirna_safe_cells(value, "identifier"), "quoted TSV")
  })
}
for (value in c("\u00a0=payload", "\u3000=payload", "MIR_A\u00a0")) {
  test("Unicode boundary whitespace is rejected without rewriting identifiers", {
    expect_error(mirna_safe_cells(value, "identifier"), "boundary whitespace")
  })
}
test("bundled tables retain exact metadata sample order", {
  fixture <- fixture_inputs()
  fixture$expression <- fixture$expression[, c(1, 7, 5, 3, 6, 2, 4)]
  inputs <- mirna_inputs(fixture$expression, fixture$samples, fixture_config())
  expect_equal(colnames(inputs$expression), fixture$samples$sample_id)
  expect_equal(nrow(inputs$expression), 30L)
  expect_equal(unname(inputs$design[, "reference"]), c(1L, 1L, 1L, 0L, 0L, 0L))
  expect_equal(inputs$collapsed_rows, 0L)
})
test("group names are data, never contrast expressions", {
  fixture <- fixture_inputs()
  fixture$samples$group <- c(rep("Reference group", 3), rep("Treatment & condition", 3))
  inputs <- mirna_inputs(fixture$expression, fixture$samples,
                          mirna_config("Reference group", "Treatment & condition", "log2"))
  expect_equal(colnames(inputs$design), c("reference", "comparison"))
})
test("synthetic replicates retain different feature rankings before quantile normalization", {
  fixture <- fixture_inputs()
  matrix <- mirna_inputs(fixture$expression, fixture$samples, fixture_config())$expression
  expect(!identical(order(matrix[, 1]), order(matrix[, 2])))
  expect(!identical(order(matrix[, 4]), order(matrix[, 5])))
})
test("duplicate features require explicit arithmetic averaging", {
  fixture <- fixture_inputs()
  expression <- rbind(fixture$expression, fixture$expression[1, ])
  expression[nrow(expression), "C1"] <- as.character(as.numeric(expression[1, "C1"]) + 2)
  expect_error(mirna_inputs(expression, fixture$samples, fixture_config()), "Duplicate miRNA")
  inputs <- mirna_inputs(expression, fixture$samples, fixture_config(duplicates = "mean"))
  expect_equal(nrow(inputs$expression), 30L)
  expect_equal(inputs$collapsed_rows, 1L)
  expect_equal(unname(inputs$expression[1, "C1"]), as.numeric(expression[1, "C1"]) + 1)
})
test("duplicate sample columns and sample metadata rows are rejected", {
  fixture <- fixture_inputs()
  bad <- fixture$expression
  names(bad)[3] <- names(bad)[2]
  expect_error(mirna_inputs(bad, fixture$samples, fixture_config()), "headers")
  fixture$samples$sample_id[2] <- fixture$samples$sample_id[1]
  expect_error(mirna_inputs(fixture$expression, fixture$samples, fixture_config()), "unique")
})
test("sample identifiers match exactly rather than by case or position", {
  fixture <- fixture_inputs()
  fixture$samples$sample_id[1] <- "c1"
  expect_error(mirna_inputs(fixture$expression, fixture$samples, fixture_config()), "match exactly")
})
test("unrepresented or extra groups and single-replicate groups are rejected", {
  fixture <- fixture_inputs()
  fixture$samples$group[1] <- "Other"
  expect_error(mirna_inputs(fixture$expression, fixture$samples, fixture_config()), "exactly the two")
  fixture$samples$group <- c("Control", rep("Treatment", 5))
  expect_error(mirna_inputs(fixture$expression, fixture$samples, fixture_config()), "at least two")
})
test("unexpected metadata covariates are not silently ignored", {
  fixture <- fixture_inputs()
  fixture$samples$subject <- letters[1:6]
  expect_error(mirna_inputs(fixture$expression, fixture$samples, fixture_config()), "exactly sample_id and group")
})
for (value in c("", "missing", "NA", "Inf", "-Inf", "NaN")) {
  test(paste("reject nonfinite expression", value), {
    fixture <- fixture_inputs()
    fixture$expression[1, 2] <- value
    expect_error(mirna_inputs(fixture$expression, fixture$samples, fixture_config()), "finite numeric")
  })
}
test("too few unique features are rejected", {
  fixture <- fixture_inputs()
  expect_error(mirna_inputs(fixture$expression[1, ], fixture$samples, fixture_config()), "At least two")
})
test("UTF-8 BOM and literal leading-zero identifiers are preserved", {
  mirna_temp(function(directory) {
    file <- file.path(directory, "input.tsv")
    writeLines(c("\ufeffmiRNA_id\tS1", "001\t1"), file, useBytes = TRUE)
    table <- mirna_read_tsv(file, "example")
    expect_equal(names(table), c("miRNA_id", "S1"))
    expect_equal(table$miRNA_id, "001")
  })
})
for (lines in list(c("A\tA", "1\t2"), c("A\tB", "1\t2\t3"),
                    c("A\tB", "1\t2", ""), c("A\tB"), c("A\tB", "1"))) {
  test("malformed TSV shape, duplicate headers and blank rows fail", {
    mirna_temp(function(directory) {
      file <- file.path(directory, "invalid.tsv")
      writeLines(lines, file)
      error <- tryCatch({mirna_read_tsv(file, "example"); NULL}, error = identity)
      expect(inherits(error, "error"))
    })
  })
}
test("missing and directory-valued inputs fail", {
  mirna_temp(function(directory) {
    expect_error(mirna_read_tsv(file.path(directory, "missing.tsv"), "example"), "Missing")
    expect_error(mirna_read_tsv(directory, "example"), "Missing")
  })
})
test("selection thresholds are strict, not inclusive", {
  raw <- fixture_statistics()
  result <- mirna_statistics(raw, rownames(raw), fixture_config())
  expect_equal(result$miRNA_id[result$selected], c("MIR_EXAMPLE_01", "MIR_EXAMPLE_02"))
  expect_equal(result$direction, c("up", "down", "not_selected", "not_selected", "not_selected"))
})
test("feature coverage, duplicate rows and probabilities are validated", {
  raw <- fixture_statistics()
  expect_error(mirna_statistics(raw[-1, ], rownames(raw), fixture_config()), "complete")
  bad <- raw
  bad$adj.P.Val[1] <- -1
  expect_error(mirna_statistics(bad, rownames(raw), fixture_config()), "nonfinite")
  bad <- raw
  bad$logFC[1] <- Inf
  expect_error(mirna_statistics(bad, rownames(raw), fixture_config()), "nonfinite")
})
test("ranking follows p values, and empty selections remain a complete table", {
  raw <- fixture_statistics()[5:1, ]
  result <- mirna_statistics(raw, rownames(raw), fixture_config(logfc = 100))
  expect_equal(result$miRNA_id, sprintf("MIR_EXAMPLE_%02d", 1:5))
  expect_equal(sum(result$selected), 0L)
})
test("reviewed targets require provenance and explicit review decisions", {
  table <- mirna_read_tsv(file.path(root, "examples", "reviewed_targets.tsv"), "targets")
  ids <- fixture_inputs()$expression$miRNA_id
  mirna_targets(table, ids)
  bad <- table
  bad$review[1] <- "yes"
  expect_error(mirna_targets(bad, ids), "include or exclude")
  bad <- table
  bad$miRNA_id[1] <- "mir_example_01"
  expect_error(mirna_targets(bad, ids), "present")
  bad <- table
  bad$reference[1] <- ""
  expect_error(mirna_targets(bad, ids), "nonblank")
})
test("target evidence preserves rows while gene lists alone are deduplicated", {
  raw <- fixture_statistics()
  statistics <- mirna_statistics(raw, rownames(raw), fixture_config())
  table <- data.frame(miRNA_id = c("MIR_EXAMPLE_01", "MIR_EXAMPLE_02", "MIR_EXAMPLE_01", "MIR_EXAMPLE_03", "MIR_EXAMPLE_01"),
                       target_gene_id = c("G_SHARED", "G_SHARED", "G_SHARED", "G_OTHER", "G_EXCLUDED"),
                       database = "synthetic", evidence_type = "illustrative", reference = "example-only",
                       review = c("include", "include", "include", "include", "exclude"))
  result <- mirna_target_exports(mirna_targets(table, statistics$miRNA_id), statistics)
  expect_equal(nrow(result$review), 5L)
  expect_equal(nrow(result$evidence), 3L)
  expect_equal(result$all, "G_SHARED")
  expect_equal(result$up, "G_SHARED")
  expect_equal(result$down, "G_SHARED")
  expect_equal(result$review$reason[4:5], c("miRNA_not_selected", "excluded_by_review"))
  expect(is.null(mirna_target_exports(NULL, statistics)))
})
test("target schema and absent statistic keys are rejected", {
  table <- mirna_read_tsv(file.path(root, "examples", "reviewed_targets.tsv"), "targets")
  expect_error(mirna_targets(table[, -1], character()), "Target evidence needs")
  raw <- fixture_statistics()
  expect_error(mirna_target_exports(table, mirna_statistics(raw, rownames(raw), fixture_config())), "does not match")
})
test("heatmaps use metadata order and exclude constant selected rows", {
  raw <- fixture_statistics()
  statistics <- mirna_statistics(raw, rownames(raw), fixture_config())
  matrix <- matrix(seq_len(30), nrow = 5, dimnames = list(rownames(raw), paste0("S", 1:6)))
  matrix[2, ] <- 2
  display <- mirna_heatmap_matrix(matrix, statistics, 40)
  expect_equal(rownames(display), "MIR_EXAMPLE_01")
  expect_equal(colnames(display), colnames(matrix))
  expect(abs(rowMeans(display)) < 1e-12)
  expect_equal(unname(apply(display, 1, sd)), 1)
})
test("plots handle zero probabilities, a single heatmap row and empty selections", {
  mirna_temp(function(directory) {
    raw <- fixture_statistics()
    raw$P.Value[1] <- 0
    raw$adj.P.Val[1] <- 0
    statistics <- mirna_statistics(raw, rownames(raw), fixture_config())
    matrix <- matrix(seq_len(30), nrow = 5, dimnames = list(rownames(raw), paste0("S", 1:6)))
    matrix[2, ] <- 2
    status <- mirna_plots(matrix, statistics, fixture_config(), directory, TRUE)
    expect(all(status$status == "completed"))
    expect(file.info(file.path(directory, "volcano.pdf"))$size > 100)
    expect(file.info(file.path(directory, "heatmap.pdf"))$size > 100)
    empty <- mirna_statistics(raw, rownames(raw), fixture_config(logfc = 100))
    status <- mirna_plots(matrix, empty, fixture_config(), directory, TRUE)
    expect_equal(status$status[2], "skipped_no_selected_variable_rows")
    expect(all(mirna_plots(matrix, empty, fixture_config(), directory, FALSE)$status == "not_requested"))
  })
})
test("existing output files/directories and missing parents are never accepted", {
  mirna_temp(function(directory) {
    expect_error(mirna_new_output(directory), "already exists")
    file <- file.path(directory, "keep.txt")
    writeLines("preserve", file)
    expect_error(mirna_new_output(file), "already exists")
    expect_error(mirna_new_output(file.path(directory, "missing", "output")), "parent")
    expect_equal(readLines(file), "preserve")
    mirna_new_output(file.path(directory, "new"))
  })
})
test("manifests omit input locations and preserve only role, size and checksum", {
  file <- file.path(root, "examples", "expression.tsv")
  manifest <- mirna_manifest(file, "expression")
  expect_equal(names(manifest), c("role", "bytes", "md5"))
  expect_equal(manifest$md5, unname(tools::md5sum(file)))
  expect(!any(grepl(root, unlist(manifest), fixed = TRUE)))
  expect(!any(grepl("\\\\", mirna_source_manifest(root)$role)))
})
test("CLI help, example defaults and explicit user arguments are consistent", {
  expect_equal(mirna_parse("--help")$action, "help")
  expect_equal(mirna_parse("--check-environment")$action, "check-environment")
  options <- mirna_parse(c("--example", "--output", "new", "--plots"))
  expect(options$example && options$plots)
  expect_equal(options$config$input_scale, "log2")
  options <- mirna_parse(c("--expression", "E", "--samples", "S", "--reference", "C",
                           "--comparison", "T", "--input-scale", "log2", "--output", "O",
                           "--normalization", "none", "--duplicates", "mean"))
  expect_equal(options$config$normalization, "none")
  expect_equal(options$config$duplicates, "mean")
})
for (arguments in list(character(), c("--unknown"), c("--output"), c("--help", "--plots"),
                       c("--check-environment", "--plots"), c("--example", "--example", "--output", "O"),
                       c("--example", "--output", "O", "--expression", "E"),
                       c("--example", "--output", "O", "--fdr", "no"),
                       c("--expression", "E", "--samples", "S", "--reference", "C", "--comparison", "T", "--output", "O"))) {
  test("invalid CLI combinations fail without analysis", {
    error <- tryCatch({mirna_parse(arguments); NULL}, error = identity)
    expect(inherits(error, "error"))
  })
}
test("missing limma has no installer, network access or substitute statistics", {
  source <- paste(readLines(file.path(root, "src", "workflow.R")), collapse = "\n")
  expect(!grepl("install.packages|BiocManager|download.file|url\\(", source))
  if (!requireNamespace("limma", quietly = TRUE)) {
    mirna_temp(function(directory) {
      output <- file.path(directory, "new")
      expect_error(mirna_run(file.path(root, "examples", "expression.tsv"),
                              file.path(root, "examples", "samples.tsv"), output,
                              fixture_config(), root), "limma is unavailable")
      expect(!dir.exists(output))
    })
  }
})
