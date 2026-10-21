array_temp <- function(code) {
  directory <- tempfile("array-check-")
  if (!dir.create(directory)) stop("Cannot create an isolated temporary directory.")
  on.exit({
    if (dir.exists(directory)) {
      target <- normalizePath(directory, winslash = "/", mustWork = TRUE)
      parent <- normalizePath(tempdir(), winslash = "/", mustWork = TRUE)
      if (!startsWith(target, paste0(parent, "/")) || !startsWith(basename(target), "array-check-")) {
        stop("Refusing cleanup outside the isolated temporary directory.")
      }
      unlink(target, recursive = TRUE)
    }
  })
  code(directory)
}

array_fixture <- function(layer = "mrna") {
  files <- array_example_files(root)
  list(expression = array_read(files[[paste0(layer, "_expression")]], "expression"),
       annotation = array_read(files[[paste0(layer, "_annotation")]], "annotation"),
       samples = array_read(files[[paste0(layer, "_samples")]], "samples"))
}

array_prepared <- function(layer = "mrna", scale = "intensity") {
  input <- array_fixture(layer)
  array_prepare(input$expression, input$annotation, input$samples, scale)
}

array_none_config <- function(plots = FALSE) array_config("intensity", "intensity", "none", "none", plots)
