test("quantile normalization agrees with an independent no-tie reference", {
  matrix <- matrix(c(5, 1, 3, 4, 2, 6, 9, 7, 8), nrow = 3,
                     dimnames = list(paste0("G", 1:3), paste0("S", 1:3)))
  expected <- matrix(c(20/3, 10/3, 5, 5, 10/3, 20/3, 20/3, 10/3, 5),
                       nrow = 3, dimnames = dimnames(matrix))
  expect_equal(array_normalize(matrix, "quantile"), expected)
})
test("tied expression values receive equal pooled-reference values", {
  matrix <- matrix(c(1, 1, 3, 1, 3, 3), nrow = 3,
                     dimnames = list(paste0("G", 1:3), paste0("S", 1:2)))
  expected <- matrix(c(1.5, 1.5, 3, 1, 2.5, 2.5), nrow = 3, dimnames = dimnames(matrix))
  expect_equal(array_normalize(matrix, "quantile"), expected)
})
test("default independent branches match direct limma normalization", {
  array_temp(function(directory) {
    files <- array_example_files(root)
    before <- tools::md5sum(files)
    output <- file.path(directory, "result")
    result <- array_run(files, output, array_config("intensity", "intensity", plots = TRUE), root)
    mrna <- array_prepared("mrna")$collapsed
    mirna <- array_prepared("mirna")$collapsed
    expected <- limma::normalizeQuantiles(mirna, ties = TRUE)
    dimnames(expected) <- dimnames(mirna)
    expect_equal(result$mrna$processed, mrna)
    expect_equal(result$mirna$processed, expected)
    expect_equal(tools::md5sum(files), before)
    written <- read.delim(file.path(output, "mirna_processed_expression.tsv"), check.names = FALSE)
    expect_equal(written$feature_id, rownames(expected))
    expect_equal(names(written)[-1], colnames(expected))
    expect_equal(unname(as.matrix(written[, -1])), unname(expected))
    metadata <- read.delim(file.path(output, "run_metadata.tsv"), check.names = FALSE)
    expect_equal(metadata$value[metadata$key == "status"], "completed")
    expect_equal(metadata$value[metadata$key == "mrna_normalization"], "none")
    expect_equal(metadata$value[metadata$key == "mirna_normalization"], "quantile")
  })
})
test("changing one assay cannot alter the other assay's processed matrix", {
  array_temp(function(directory) {
    files <- array_example_files(root)
    baseline <- array_run(files, file.path(directory, "first"),
                            array_config("intensity", "intensity"), root)
    changed <- array_read(files[["mirna_expression"]], "expression")
    changed[, -1] <- lapply(changed[, -1], function(values) as.character(as.numeric(values) * 10 + 4))
    altered <- file.path(directory, "altered.tsv")
    array_write(changed, altered)
    files[["mirna_expression"]] <- altered
    second <- array_run(files, file.path(directory, "second"),
                          array_config("intensity", "intensity"), root)
    expect_equal(second$mrna$processed, baseline$mrna$processed)
    expect(!isTRUE(all.equal(second$mirna$processed, baseline$mirna$processed)))
  })
})
