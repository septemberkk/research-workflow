mirna_temp <- function(code) {
  directory <- tempfile("mirna-check-")
  if (!dir.create(directory)) stop("Cannot create a temporary check directory.")
  on.exit({
    if (dir.exists(directory)) {
      resolved <- normalizePath(directory, winslash = "/", mustWork = TRUE)
      parent <- normalizePath(tempdir(), winslash = "/", mustWork = TRUE)
      if (!startsWith(resolved, paste0(parent, "/")) || !startsWith(basename(resolved), "mirna-check-")) {
        stop("Refusing to clean an unexpected temporary directory.")
      }
      unlink(resolved, recursive = TRUE)
    }
  })
  code(directory)
}

fixture_config <- function(...) mirna_config("Control", "Treatment", "log2", ...)

fixture_statistics <- function() {
  data.frame(logFC = c(2, -2, 1, 3, 0), AveExpr = c(5, 5, 5, 5, 5),
              t = c(4, -4, 3, 2, 0), P.Value = c(.001, .002, .003, .01, 1),
              adj.P.Val = c(.005, .01, .02, .05, 1), B = c(2, 1, 0, -1, -5),
              row.names = sprintf("MIR_EXAMPLE_%02d", 1:5))
}

fixture_inputs <- function() {
  list(expression = mirna_read_tsv(file.path(root, "examples", "expression.tsv"), "expression"),
       samples = mirna_read_tsv(file.path(root, "examples", "samples.tsv"), "samples"))
}
