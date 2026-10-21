array_fail <- function(message) stop(message, call. = FALSE)

array_text <- function(value, label, empty = FALSE) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      (!empty && !nzchar(value)) ||
      value != trimws(value, whitespace = "[\\h\\v]") || grepl("[[:cntrl:]]", value)) {
    array_fail(paste0(label, " must be plain single-line text without boundary whitespace."))
  }
  value
}

array_cells <- function(values, label, empty = FALSE) {
  for (value in values) {
    array_text(value, label, empty)
    if (grepl('"', value, fixed = TRUE) || grepl("^[=+@-]", value) ||
        grepl("\\p{Cf}|\\p{Zl}|\\p{Zp}", enc2utf8(value), perl = TRUE)) {
      array_fail(paste0(label, " contains unsafe quoting, formula-like or invisible text."))
    }
  }
  invisible(values)
}

array_read <- function(file, label) {
  array_text(file, label)
  if (!file.exists(file) || dir.exists(file)) array_fail(paste0("Missing ", label, " input."))
  lines <- tryCatch(readLines(file, encoding = "UTF-8", warn = FALSE),
                    error = function(error) array_fail(paste0("Cannot read ", label, " input.")))
  if (length(lines) < 2L || any(!validUTF8(lines)) || any(!nzchar(lines))) {
    array_fail("Inputs must be UTF-8 TSV with data rows and no blank physical lines.")
  }
  lines[1L] <- sub("^\ufeff", "", lines[1L])
  con <- textConnection(lines)
  on.exit(close(con))
  widths <- count.fields(con, sep = "\t", quote = "", comment.char = "", blank.lines.skip = FALSE)
  if (length(widths) != length(lines) || widths[1L] < 2L || any(widths != widths[1L])) {
    array_fail("Inconsistent tab-separated field counts.")
  }
  table <- tryCatch(read.table(text = lines, header = TRUE, sep = "\t", quote = "",
                                comment.char = "", colClasses = "character", check.names = FALSE,
                                na.strings = NULL, stringsAsFactors = FALSE, fill = FALSE),
                     error = function(error) array_fail("Invalid TSV input."))
  array_cells(names(table), "headers")
  if (!nrow(table) || anyDuplicated(names(table))) array_fail("Headers must be unique and data rows nonempty.")
  table
}

array_write <- function(table, file) {
  write.table(table, file, sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = TRUE, na = "", fileEncoding = "UTF-8")
}

array_manifest <- function(files, roles) {
  data.frame(role = roles, bytes = unname(file.info(files)$size),
              md5 = unname(tools::md5sum(files)), stringsAsFactors = FALSE)
}

array_output <- function(output) {
  array_text(output, "output")
  if (file.exists(output) || dir.exists(output)) array_fail("Output already exists; no files will be overwritten.")
  if (!dir.exists(dirname(output))) array_fail("Output parent must already exist.")
}

array_config <- function(mrna_scale, mirna_scale, mrna_normalization = "none",
                          mirna_normalization = "quantile", plots = FALSE) {
  for (value in list(mrna_scale, mirna_scale)) {
    if (length(value) != 1L || !value %in% c("intensity", "log2")) {
      array_fail("Every layer needs an explicit intensity or log2 scale.")
    }
  }
  for (value in list(mrna_normalization, mirna_normalization)) {
    if (length(value) != 1L || !value %in% c("none", "quantile")) {
      array_fail("Normalization must be none or quantile.")
    }
  }
  if (!is.logical(plots) || length(plots) != 1L || is.na(plots)) array_fail("Plots must be TRUE or FALSE.")
  list(mrna = list(scale = mrna_scale, normalization = mrna_normalization),
       mirna = list(scale = mirna_scale, normalization = mirna_normalization), plots = plots)
}
