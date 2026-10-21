array_prepare <- function(expression, annotation, samples, scale) {
  if (!is.data.frame(expression) || ncol(expression) < 3L || names(expression)[1L] != "probe_id") {
    array_fail("Expression requires probe_id and at least two sample columns.")
  }
  if (!is.data.frame(annotation) || !identical(names(annotation), c("probe_id", "feature_id"))) {
    array_fail("Annotation requires exactly probe_id and feature_id.")
  }
  if (!is.data.frame(samples) || !identical(names(samples), c("sample_id", "group", "include"))) {
    array_fail("Samples require exactly sample_id, group and include.")
  }
  array_cells(names(expression), "expression headers")
  array_cells(expression$probe_id, "expression probe IDs")
  array_cells(annotation$probe_id, "annotation probe IDs")
  array_cells(annotation$feature_id, "annotation feature IDs", empty = TRUE)
  for (column in names(samples)) array_cells(samples[[column]], paste0("sample ", column))
  if (anyDuplicated(names(expression)) || anyDuplicated(expression$probe_id) ||
      anyDuplicated(annotation$probe_id) || anyDuplicated(samples$sample_id)) {
    array_fail("Expression probes, annotation probes, sample IDs and sample columns must each be unique.")
  }
  if (!setequal(names(expression)[-1L], samples$sample_id)) {
    array_fail("Expression columns and sample metadata must match exactly.")
  }
  if (any(!samples$include %in% c("yes", "no"))) array_fail("Sample include values must be yes or no.")
  selected <- samples[samples$include == "yes", , drop = FALSE]
  if (nrow(selected) < 2L) array_fail("At least two samples must be explicitly included in each layer.")
  values <- suppressWarnings(matrix(as.numeric(as.matrix(expression[, -1L, drop = FALSE])),
                                      nrow = nrow(expression),
                                      dimnames = list(expression$probe_id, names(expression)[-1L])))
  if (any(!is.finite(values))) array_fail("Expression must contain finite numeric values; missing data are not replaced by zero.")
  if (!scale %in% c("intensity", "log2")) array_fail("Unsupported expression scale.")
  if (scale == "intensity" && any(values < 0)) array_fail("Intensity-scale values cannot be negative.")
  index <- match(expression$probe_id, annotation$probe_id)
  features <- annotation$feature_id[index]
  features[is.na(features)] <- ""
  status <- ifelse(nzchar(features), "mapped", "unmapped")
  ambiguous <- grepl("///|;|\\||[\\h\\v]", enc2utf8(features), perl = TRUE)
  status[ambiguous] <- "ambiguous"
  mapping <- data.frame(probe_id = expression$probe_id, feature_id = features, status = status)
  retained <- status == "mapped"
  if (!any(retained)) array_fail("No unambiguous mapped probes remain.")
  ids <- sort(unique(features[retained]), method = "radix")
  values <- values[, selected$sample_id, drop = FALSE]
  collapsed <- do.call(rbind, lapply(ids, function(id) {
    colMeans(values[retained & features == id, , drop = FALSE])
  }))
  dimnames(collapsed) <- list(ids, selected$sample_id)
  if (any(!is.finite(collapsed))) array_fail("Arithmetic averaging produced nonfinite values.")
  members <- data.frame(feature_id = features[retained], probe_id = expression$probe_id[retained])
  counts <- data.frame(feature_id = ids, n_probes = tabulate(match(features[retained], ids), nbins = length(ids)))
  summary <- c(input_probes = nrow(expression), mapped_probes = sum(retained),
                unmapped_probes = sum(status == "unmapped"), ambiguous_probes = sum(status == "ambiguous"),
                unused_annotation_probes = sum(!annotation$probe_id %in% expression$probe_id),
                feature_rows = length(ids), input_samples = nrow(samples),
                selected_samples = nrow(selected), excluded_samples = sum(samples$include == "no"))
  list(collapsed = collapsed, mapping = mapping, members = members, counts = counts,
       samples = samples, selected = selected, annotation = annotation,
       summary = data.frame(metric = names(summary), value = unname(summary)))
}

array_require_limma <- function() {
  if (!requireNamespace("limma", quietly = TRUE)) {
    array_fail("limma is required for quantile normalization. No packages were installed.")
  }
}

array_normalize <- function(matrix, method) {
  if (identical(method, "none")) return(matrix)
  if (!identical(method, "quantile")) array_fail("Unsupported normalization method.")
  if (nrow(matrix) < 2L) array_fail("Quantile normalization requires at least two retained feature rows.")
  array_require_limma()
  result <- limma::normalizeQuantiles(matrix, ties = TRUE)
  dimnames(result) <- dimnames(matrix)
  if (any(!is.finite(result))) array_fail("Quantile normalization produced nonfinite values.")
  result
}

array_distribution <- function(matrix, stage) {
  summary <- t(vapply(seq_len(ncol(matrix)), function(index) {
    values <- matrix[, index]
    c(minimum = min(values), q1 = unname(quantile(values, .25)),
      median = median(values), mean = mean(values), q3 = unname(quantile(values, .75)), maximum = max(values))
  }, numeric(6)))
  data.frame(sample_id = colnames(matrix), stage = stage, summary, row.names = NULL, check.names = FALSE)
}

array_plot <- function(layer, result, scale, file) {
  grDevices::pdf(file, width = 9, height = 7)
  on.exit(grDevices::dev.off())
  graphics::par(mfrow = c(2, 1), mar = c(5, 4, 3, 1))
  groups <- result$selected$group
  palette <- grDevices::hcl.colors(length(unique(groups)), "Dark 3")
  colors <- palette[match(groups, unique(groups))]
  for (stage in c("collapsed", "processed")) {
    graphics::boxplot(as.data.frame(result[[stage]], check.names = FALSE),
                       col = colors, outline = FALSE, las = 2,
                       ylab = paste("Expression:", scale), main = paste(layer, stage))
  }
}
