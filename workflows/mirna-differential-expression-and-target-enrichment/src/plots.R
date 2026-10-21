mirna_heatmap_matrix <- function(expression, statistics, limit) {
  ids <- head(statistics$miRNA_id[statistics$selected], limit)
  selected <- expression[ids, , drop = FALSE]
  if (!nrow(selected)) return(selected)
  variability <- apply(selected, 1L, stats::sd)
  selected <- selected[is.finite(variability) & variability > 0, , drop = FALSE]
  if (!nrow(selected)) return(selected)
  t(scale(t(selected), center = TRUE, scale = TRUE))
}

mirna_volcano <- function(statistics, config, file) {
  grDevices::pdf(file, width = 6, height = 5)
  on.exit(grDevices::dev.off())
  y <- -log10(pmax(statistics$adj.P.Val, .Machine$double.xmin))
  colors <- c(not_selected = "#A9A5B2", up = "#A64D79", down = "#5675A6")[statistics$direction]
  graphics::plot(statistics$logFC, y, pch = 16, cex = .6, col = colors,
                 xlab = "log2 fold change (comparison - reference)", ylab = "-log10(BH FDR)")
  graphics::abline(v = c(-config$logfc, config$logfc),
                   h = -log10(config$fdr), lty = 2, col = "#747078")
}

mirna_heatmap <- function(matrix, file) {
  if (!nrow(matrix)) return(invisible(NULL))
  grDevices::pdf(file, width = 8, height = max(4, min(12, nrow(matrix) * .2 + 2)))
  on.exit(grDevices::dev.off())
  colors <- grDevices::colorRampPalette(c("#4E6F9D", "#FAF7F2", "#A64D79"))(101)
  bound <- max(abs(matrix))
  graphics::image(seq_len(ncol(matrix)), seq_len(nrow(matrix)),
                  t(matrix[nrow(matrix):1L, , drop = FALSE]), col = colors,
                  zlim = c(-bound, bound), axes = FALSE, xlab = "", ylab = "")
  graphics::axis(1, at = seq_len(ncol(matrix)), labels = colnames(matrix), las = 2, cex.axis = .7)
  graphics::axis(2, at = seq_len(nrow(matrix)), labels = rev(rownames(matrix)), las = 2, cex.axis = .6)
  graphics::title(main = "Selected miRNAs: row z-scores")
}

mirna_plots <- function(expression, statistics, config, output, requested) {
  status <- data.frame(plot = c("volcano", "heatmap"),
                       status = c("not_requested", "not_requested"))
  if (!requested) return(status)
  mirna_volcano(statistics, config, file.path(output, "volcano.pdf"))
  status$status[1L] <- "completed"
  matrix <- mirna_heatmap_matrix(expression, statistics, config$heatmap_rows)
  if (!nrow(matrix)) {
    status$status[2L] <- "skipped_no_selected_variable_rows"
  } else {
    mirna_heatmap(matrix, file.path(output, "heatmap.pdf"))
    status$status[2L] <- "completed"
  }
  status
}
