mirna_targets <- function(table, assayed_ids) {
  required <- c("miRNA_id", "target_gene_id", "database", "evidence_type", "reference", "review")
  if (!is.data.frame(table) || !identical(names(table), required)) {
    mirna_fail("Target evidence needs miRNA_id, target_gene_id, database, evidence_type, reference and review.")
  }
  for (column in names(table)) mirna_safe_cells(table[[column]], paste0("target ", column))
  if (any(!table$review %in% c("include", "exclude"))) {
    mirna_fail("Target review values must be include or exclude.")
  }
  if (any(!table$miRNA_id %in% assayed_ids)) {
    mirna_fail("Target miRNA identifiers must be present in the assayed expression matrix.")
  }
  table
}

mirna_target_exports <- function(targets, statistics) {
  if (is.null(targets)) return(NULL)
  index <- match(targets$miRNA_id, statistics$miRNA_id)
  if (anyNA(index)) mirna_fail("Target evidence does not match the statistics.")
  result <- targets
  result$miRNA_logFC <- statistics$logFC[index]
  result$miRNA_BH_FDR <- statistics$adj.P.Val[index]
  result$miRNA_direction <- statistics$direction[index]
  result$selected_for_followup <- result$review == "include" & statistics$selected[index]
  result$reason <- ifelse(result$review == "exclude", "excluded_by_review",
                          ifelse(statistics$selected[index], "selected", "miRNA_not_selected"))
  included <- result[result$selected_for_followup, , drop = FALSE]
  gene_list <- function(direction = NULL) {
    rows <- if (is.null(direction)) included else included[included$miRNA_direction == direction, , drop = FALSE]
    sort(unique(rows$target_gene_id), method = "radix")
  }
  list(review = result, evidence = included, all = gene_list(),
       up = gene_list("up"), down = gene_list("down"))
}
