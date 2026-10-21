test("branch settings preserve explicit scales and independent defaults", {
  config <- array_config("intensity", "log2")
  expect_equal(config$mrna$normalization, "none")
  expect_equal(config$mirna$normalization, "quantile")
  expect_equal(config$mirna$scale, "log2")
})
for (value in list("auto", "counts", "", NULL, NA_character_)) {
  test("undeclared or unsupported scales are rejected", {
    expect_error(array_config(value, "intensity"), "explicit")
  })
}
for (value in list("auto", "", NULL, NA_character_)) {
  test("unsupported normalization is rejected", {
    expect_error(array_config("intensity", "intensity", value), "Normalization")
  })
}
test("plot request must be scalar logical", {
  expect_error(array_none_config(plots = "yes"), "Plots")
  expect_error(array_none_config(plots = NA), "Plots")
})
for (value in c(" A", "A ", "\u00a0A", "A\u3000", "A\nB", "A\u0085B", "")) {
  test("invalid whitespace or blank identifiers are rejected without trimming", {
    expect_error(array_cells(value, "identifier"), "plain")
  })
}
for (value in c('"=payload"', "=payload", "+payload", "@payload", "-payload", "A\u200b", "A\u00ad", "A\u061c", "A\u2028B", "A\u2029B")) {
  test("quoted, formula-like and Unicode formatting text cannot enter exports", {
    expect_error(array_cells(value, "identifier"), "unsafe")
  })
}
test("numeric-looking and missing-value-looking strings remain literal IDs", {
  array_cells(c("001", "NA", "MIR_EXAMPLE", "GENE_\u03b1"), "identifier")
})
test("selected columns follow explicit metadata order", {
  result <- array_prepared()
  expect_equal(colnames(result$collapsed), c("M_C2", "M_T2", "M_C1", "M_T1"))
  expect_equal(result$selected$sample_id, colnames(result$collapsed))
  expect_equal(ncol(result$collapsed), 4L)
})
test("mapped probes are averaged equally on the supplied scale", {
  result <- array_prepared()
  expect_equal(unname(result$collapsed["GENE_EXAMPLE_01", ]), c(14, 23, 12, 22))
  expect_equal(unname(result$collapsed["GENE_EXAMPLE_02", ]), c(5, 9, 4, 6))
  expect_equal(result$counts$n_probes, c(2L, 2L, 2L))
  expect_equal(nrow(result$members), 6L)
})
test("every expression probe retains an explicit mapping outcome", {
  result <- array_prepared()
  expect_equal(nrow(result$mapping), 8L)
  expect_equal(result$mapping$status, c(rep("mapped", 6), "ambiguous", "unmapped"))
  expect_equal(result$summary$value[result$summary$metric == "ambiguous_probes"], 1L)
  expect_equal(result$summary$value[result$summary$metric == "excluded_samples"], 2L)
})
test("a blank annotation stays unmapped, never a guessed feature", {
  fixture <- array_fixture()
  fixture$annotation$feature_id[1] <- ""
  result <- array_prepare(fixture$expression, fixture$annotation, fixture$samples, "intensity")
  expect_equal(result$mapping$status[1], "unmapped")
  expect_equal(result$mapping$feature_id[1], "")
})
for (value in c("A///B", "A;B", "A|B", "A B", "A\u00a0B")) {
  test("compound mappings are reported instead of taking their first symbol", {
    fixture <- array_fixture()
    fixture$annotation$feature_id[1] <- value
    result <- array_prepare(fixture$expression, fixture$annotation, fixture$samples, "intensity")
    expect_equal(result$mapping$status[1], "ambiguous")
    expect_equal(result$mapping$feature_id[1], value)
  })
}
test("unused annotation entries remain separate from observed probes", {
  fixture <- array_fixture()
  fixture$annotation <- rbind(fixture$annotation, data.frame(probe_id = "UNOBSERVED", feature_id = "GENE_OTHER"))
  result <- array_prepare(fixture$expression, fixture$annotation, fixture$samples, "intensity")
  expect_equal(result$summary$value[result$summary$metric == "unused_annotation_probes"], 1L)
  expect(!"GENE_OTHER" %in% rownames(result$collapsed))
})
test("expression-column order does not replace the metadata order", {
  fixture <- array_fixture()
  fixture$expression <- fixture$expression[, c(1, 7, 6, 5, 4, 3, 2)]
  result <- array_prepare(fixture$expression, fixture$annotation, fixture$samples, "intensity")
  expect_equal(result$collapsed, array_prepared()$collapsed)
})
test("excluded sample values cannot affect selected feature averages", {
  fixture <- array_fixture()
  fixture$expression$M_X1 <- "100000"
  fixture$expression$M_X2 <- "200000"
  result <- array_prepare(fixture$expression, fixture$annotation, fixture$samples, "intensity")
  expect_equal(result$collapsed, array_prepared()$collapsed)
})
test("duplicate expression probes and annotation keys are rejected", {
  fixture <- array_fixture()
  expect_error(array_prepare(rbind(fixture$expression, fixture$expression[1, ]),
                              fixture$annotation, fixture$samples, "intensity"), "unique")
  expect_error(array_prepare(fixture$expression, rbind(fixture$annotation, fixture$annotation[1, ]),
                              fixture$samples, "intensity"), "unique")
})
test("duplicate metadata, missing metadata and case-only matches are rejected", {
  fixture <- array_fixture()
  expect_error(array_prepare(fixture$expression, fixture$annotation,
                              rbind(fixture$samples, fixture$samples[1, ]), "intensity"), "unique")
  bad <- fixture$samples
  bad$sample_id[1] <- "m_c2"
  expect_error(array_prepare(fixture$expression, fixture$annotation, bad, "intensity"), "match exactly")
})
test("unknown include flags and fewer than two included samples are rejected", {
  fixture <- array_fixture()
  bad <- fixture$samples
  bad$include[1] <- "true"
  expect_error(array_prepare(fixture$expression, fixture$annotation, bad, "intensity"), "yes or no")
  bad$include <- c("yes", rep("no", 5))
  expect_error(array_prepare(fixture$expression, fixture$annotation, bad, "intensity"), "At least two")
})
test("preparation does not require a differential-expression group design", {
  fixture <- array_fixture()
  fixture$samples$group <- "OneGroup"
  result <- array_prepare(fixture$expression, fixture$annotation, fixture$samples, "intensity")
  expect_equal(unique(result$selected$group), "OneGroup")
})
for (value in c("", "missing", "NA", "Inf", "-Inf", "NaN")) {
  test("nonnumeric and nonfinite expression is rejected rather than zero-filled", {
    fixture <- array_fixture()
    fixture$expression[1, 2] <- value
    expect_error(array_prepare(fixture$expression, fixture$annotation, fixture$samples, "intensity"), "finite numeric")
  })
}
test("negative values require an explicitly declared log2 scale", {
  fixture <- array_fixture()
  fixture$expression[1, 2] <- "-3"
  expect_error(array_prepare(fixture$expression, fixture$annotation, fixture$samples, "intensity"), "cannot be negative")
  result <- array_prepare(fixture$expression, fixture$annotation, fixture$samples, "log2")
  expect_equal(unname(result$collapsed["GENE_EXAMPLE_01", "M_C1"]), 5.5)
})
test("an all-unmapped layer does not generate substitute results", {
  fixture <- array_fixture()
  fixture$annotation$feature_id <- ""
  expect_error(array_prepare(fixture$expression, fixture$annotation, fixture$samples, "intensity"), "No unambiguous")
})
test("none normalization preserves matrix values and dimnames", {
  matrix <- array_prepared()$collapsed
  expect(identical(array_normalize(matrix, "none"), matrix))
  expect_error(array_normalize(matrix, "auto"), "Unsupported")
})
test("summary statistics use the actual prepared sample values", {
  matrix <- array_prepared()$collapsed
  summary <- array_distribution(matrix, "collapsed")
  expect_equal(summary$sample_id, colnames(matrix))
  expect_equal(summary$mean, unname(colMeans(matrix)))
  expect_equal(summary$median, unname(apply(matrix, 2, median)))
})
test("UTF-8 BOM and leading-zero IDs survive TSV reading", {
  array_temp(function(directory) {
    file <- file.path(directory, "input.tsv")
    writeLines(c("\ufeffprobe_id\tS1", "001\t2"), file, useBytes = TRUE)
    result <- array_read(file, "expression")
    expect_equal(result$probe_id, "001")
    expect_equal(result$S1, "2")
  })
})
for (lines in list(c("A\tA", "1\t2"), c("A\tB", "1\t2\t3"),
                    c("A\tB", "1\t2", ""), c("A\tB"), c("A\tB", "1"))) {
  test("malformed TSV structure is rejected", {
    array_temp(function(directory) {
      file <- file.path(directory, "invalid.tsv")
      writeLines(lines, file)
      error <- tryCatch({array_read(file, "input"); NULL}, error = identity)
      expect(inherits(error, "error"))
    })
  })
}
test("missing inputs and existing outputs never cause rewriting", {
  array_temp(function(directory) {
    expect_error(array_read(file.path(directory, "missing"), "input"), "Missing")
    expect_error(array_output(directory), "already exists")
    file <- file.path(directory, "keep.txt")
    writeLines("preserve", file)
    expect_error(array_output(file), "already exists")
    expect_equal(readLines(file), "preserve")
    expect_error(array_output(file.path(directory, "absent", "result")), "parent")
  })
})
test("an end-to-end no-normalization run is independent and path-free", {
  array_temp(function(directory) {
    files <- array_example_files(root)
    before <- tools::md5sum(files)
    output <- file.path(directory, "result")
    results <- array_run(files, output, array_none_config(plots = TRUE), root)
    expect_equal(results$mrna$processed, array_prepared()$collapsed)
    expect_equal(results$mirna$processed, array_prepared("mirna")$collapsed)
    expect_equal(tools::md5sum(files), before)
    metadata <- read.delim(file.path(output, "run_metadata.tsv"), check.names = FALSE)
    expect_equal(metadata$value[metadata$key == "status"], "completed")
    for (key in c("differential_expression", "cross_layer_integration")) {
      expect_equal(metadata$value[metadata$key == key], "not_run")
    }
    expect_equal(metadata$value[metadata$key == "cross_layer_pairing"], "not_inferred")
    for (file in c("run_metadata.tsv", "input_manifest.tsv")) {
      text <- readLines(file.path(output, file))
      expect(!any(grepl(root, text, fixed = TRUE)))
      expect(!any(grepl(directory, text, fixed = TRUE)))
    }
    manifest <- read.delim(file.path(output, "output_manifest.tsv"), check.names = FALSE)
    expect_equal(unname(tools::md5sum(file.path(output, manifest$role))), manifest$md5)
    expect(file.info(file.path(output, "mrna_distributions.pdf"))$size > 100)
    expect(file.info(file.path(output, "mirna_distributions.pdf"))$size > 100)
    expect_error(array_run(files, output, array_none_config(), root), "already exists")
  })
})
test("an invalid second layer creates no partial first-layer output", {
  array_temp(function(directory) {
    files <- array_example_files(root)
    files[["mirna_samples"]] <- file.path(directory, "absent")
    output <- file.path(directory, "result")
    expect_error(array_run(files, output, array_none_config(), root), "missing")
    expect(!dir.exists(output))
  })
})
test("CLI requires explicit inputs and scales without accepting ambiguous combinations", {
  expect_equal(array_parse("--help")$action, "help")
  expect_equal(array_parse("--check-environment")$action, "check-environment")
  options <- array_parse(c("--example", "--output", "new", "--mirna-normalization", "none"))
  expect(options$example)
  expect_equal(options$config$mirna$normalization, "none")
})
for (request in list(character(), "--unknown", "--output", c("--help", "--plots"),
                   c("--example", "--example"), c("--example", "--output", "O", "--mrna-scale", "log2"),
                   c("--example", "--output", "O", "--mrna-normalization", "auto"))) {
  test("invalid CLI requests fail before processing", {
    error <- tryCatch({array_parse(request); NULL}, error = identity)
    expect(inherits(error, "error"))
  })
}
test("a standalone example CLI runs without limma when normalization is disabled", {
  array_temp(function(directory) {
    executable <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
    output <- file.path(directory, "cli")
    log <- system2(executable, c("--vanilla", shQuote(file.path(root, "run_all.R")), "--example",
                                  "--output", shQuote(output), "--mirna-normalization", "none"),
                     stdout = TRUE, stderr = TRUE)
    expect(is.null(attr(log, "status")))
    expect(file.exists(file.path(output, "mrna_processed_expression.tsv")))
    expect(file.exists(file.path(output, "mirna_processed_expression.tsv")))
  })
})
test("quantile normalization requires at least two retained features", {
  expect_error(array_normalize(matrix(c(1, 2), nrow = 1), "quantile"), "at least two")
})

test("missing limma never triggers an installer or substitute quantile algorithm", {
  if (!requireNamespace("limma", quietly = TRUE)) {
    expect_error(array_normalize(array_prepared()$collapsed, "quantile"), "limma is required")
    array_temp(function(directory) {
      output <- file.path(directory, "result")
      expect_error(array_run(array_example_files(root), output,
                               array_config("intensity", "intensity"), root), "limma is required")
      expect(!dir.exists(output))
    })
  }
  source <- paste(unlist(lapply(c("io.R", "prepare.R", "workflow.R", "cli.R"), function(file) {
    readLines(file.path(root, "src", file))
  })), collapse = "\n")
  expect(!grepl("install.packages|BiocManager|download.file|https?://", source))
})
