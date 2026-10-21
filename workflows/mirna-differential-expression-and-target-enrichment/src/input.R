mirna_fail <- function(message) stop(message, call. = FALSE)

mirna_text <- function(value, label, empty = FALSE) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      (!empty && !nzchar(value)) || value != trimws(value) ||
      grepl("[[:cntrl:]]", value)) {
    mirna_fail(paste0(label, " must be nonblank single-line text without boundary whitespace."))
  }
  value
}

mirna_safe_cells <- function(values, label, empty = FALSE) {
  for (value in values) {
    mirna_text(value, label, empty = empty)
    if (grepl("^[=+@-]", value) ||
        grepl("[\u200b-\u200f\u202a-\u202e\u2060-\u206f\ufeff]", value, perl = TRUE)) {
      mirna_fail(paste0(label, " contains formula-like or invisible formatting text."))
    }
  }
  invisible(values)
}

mirna_number <- function(value, label, minimum, maximum = Inf, integer = FALSE) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < minimum || value > maximum || (integer && value != floor(value))) {
    mirna_fail(paste0("Invalid ", label, "."))
  }
  value
}

mirna_config <- function(reference, comparison, input_scale, normalization = "quantile",
                         duplicates = "error", fdr = 0.05, logfc = 1, heatmap_rows = 40) {
  mirna_safe_cells(reference, "reference")
  mirna_safe_cells(comparison, "comparison")
  if (length(reference) != 1L || length(comparison) != 1L || identical(reference, comparison)) {
    mirna_fail("Reference and comparison must be different scalar group labels.")
  }
  if (!identical(input_scale, "log2")) mirna_fail("Only explicitly declared log2 expression is supported.")
  if (length(normalization) != 1L || !normalization %in% c("quantile", "none")) {
    mirna_fail("Normalization must be quantile or none.")
  }
  if (length(duplicates) != 1L || !duplicates %in% c("error", "mean")) {
    mirna_fail("Duplicate policy must be error or mean.")
  }
  mirna_number(fdr, "FDR threshold", .Machine$double.eps, 1)
  mirna_number(logfc, "logFC threshold", 0)
  mirna_number(heatmap_rows, "heatmap row limit", 1, 1000, integer = TRUE)
  list(reference = reference, comparison = comparison, input_scale = input_scale,
       normalization = normalization, duplicates = duplicates, fdr = fdr,
       logfc = logfc, heatmap_rows = as.integer(heatmap_rows))
}

mirna_read_tsv <- function(file, label) {
  mirna_text(file, label)
  if (!file.exists(file) || dir.exists(file)) mirna_fail(paste0("Missing ", label, " file."))
  lines <- tryCatch(readLines(file, encoding = "UTF-8", warn = FALSE),
                    error = function(error) mirna_fail(paste0("Cannot read ", label, ".")))
  if (length(lines) < 2L || any(!validUTF8(lines)) || any(!nzchar(lines))) {
    mirna_fail(paste0(label, " must be a nonempty UTF-8 table without blank physical lines."))
  }
  lines[1L] <- sub("^\ufeff", "", lines[1L])
  connection <- textConnection(lines)
  on.exit(close(connection))
  widths <- count.fields(connection, sep = "\t", quote = "", comment.char = "",
                         blank.lines.skip = FALSE)
  if (length(widths) != length(lines) || any(widths != widths[1L]) || widths[1L] < 2L) {
    mirna_fail(paste0(label, " has inconsistent tab-separated field counts."))
  }
  result <- tryCatch(
    read.table(text = lines, header = TRUE, sep = "\t", quote = "", comment.char = "",
               colClasses = "character", check.names = FALSE, na.strings = NULL,
               stringsAsFactors = FALSE, fill = FALSE, blank.lines.skip = FALSE),
    error = function(error) mirna_fail(paste0("Invalid ", label, " table."))
  )
  if (!nrow(result) || anyDuplicated(names(result))) mirna_fail(paste0("Invalid ", label, " headers or rows."))
  mirna_safe_cells(names(result), paste0(label, " headers"))
  result
}

mirna_inputs <- function(expression, samples, config) {
  if (!is.data.frame(expression) || ncol(expression) < 5L ||
      !identical(names(expression)[1L], "miRNA_id")) {
    mirna_fail("Expression needs miRNA_id followed by at least four sample columns.")
  }
  if (!is.data.frame(samples) || !identical(names(samples), c("sample_id", "group"))) {
    mirna_fail("Samples must contain exactly sample_id and group.")
  }
  mirna_safe_cells(names(expression), "expression headers")
  if (anyDuplicated(names(expression))) mirna_fail("Expression headers must be unique.")
  mirna_safe_cells(expression$miRNA_id, "miRNA identifiers")
  mirna_safe_cells(samples$sample_id, "sample identifiers")
  mirna_safe_cells(samples$group, "group labels")
  if (anyDuplicated(samples$sample_id)) mirna_fail("Sample identifiers must be unique.")
  if (!setequal(names(expression)[-1L], samples$sample_id)) {
    mirna_fail("Expression columns and sample identifiers must match exactly.")
  }
  if (!setequal(unique(samples$group), c(config$reference, config$comparison))) {
    mirna_fail("Samples must describe exactly the two specified independent groups.")
  }
  if (any(table(factor(samples$group, levels = c(config$reference, config$comparison))) < 2L)) {
    mirna_fail("Each group needs at least two independent biological samples.")
  }
  values <- expression[, samples$sample_id, drop = FALSE]
  numeric <- suppressWarnings(matrix(as.numeric(as.matrix(values)), nrow = nrow(values),
                                     dimnames = list(expression$miRNA_id, samples$sample_id)))
  if (any(!is.finite(numeric))) mirna_fail("Expression values must be finite numeric log2 values.")
  if (anyDuplicated(expression$miRNA_id)) {
    if (config$duplicates == "error") mirna_fail("Duplicate miRNA identifiers require an explicit mean policy.")
    ids <- unique(expression$miRNA_id)
    numeric <- do.call(rbind, lapply(ids, function(id) {
      colMeans(numeric[rownames(numeric) == id, , drop = FALSE])
    }))
    dimnames(numeric) <- list(ids, samples$sample_id)
  }
  if (nrow(numeric) < 2L) mirna_fail("At least two unique miRNA features are required.")
  design <- cbind(reference = as.integer(samples$group == config$reference),
                  comparison = as.integer(samples$group == config$comparison))
  rownames(design) <- samples$sample_id
  list(expression = numeric, samples = samples, design = design,
       input_rows = nrow(expression), collapsed_rows = nrow(expression) - nrow(numeric))
}

mirna_write_tsv <- function(value, file) {
  write.table(value, file, sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = TRUE, na = "", fileEncoding = "UTF-8")
}

mirna_write_ids <- function(ids, file) {
  writeLines(ids, file, useBytes = TRUE)
}

mirna_new_output <- function(output) {
  mirna_text(output, "output")
  if (file.exists(output) || dir.exists(output)) mirna_fail("Output already exists; choose a new directory.")
  if (!dir.exists(dirname(output))) mirna_fail("The output parent directory must already exist.")
  invisible(output)
}
