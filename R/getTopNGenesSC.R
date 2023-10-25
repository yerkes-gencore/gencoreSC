#' Return top genes from Seurat::FindMarkers result
#'
#' Returns top N genes from a `Seurat::Findmarkers` or other differential expression
#' result, based on some optional filtering criteria. Mainly intended for quick
#' returns of interesting genes to plot.
#'
#' @param result data.frame from `Seurat::FindMarkers` or similar function
#' @param N Number of genes to return, max
#' @param min_logFC Minimum absolute value for logFC
#' @param min_padj Minimum adjusted p value
#' @param min_pct Minimum percent expressed
#' @param exclusion_patterns Patterns to exclude. Defaults to ENS IDs, mitochondrial
#'  genes, and ribosomal genes
#' @param direction one of "mixed, up, down, equal". Dictates the direction
#'  of expression for returned genes. 'Mixed' does no filtering based on direction.
#'  Up and down specify positive and negative log FC. Equal returns an equal
#'  number of up and down regulated genes.
#' @param lfc_col Column containing log-fold change values
#' @param pval_col Column containing adjusted p values
#'
#' @returns A vector of gene names
#' @export
#'
#' @examples
#' \dontrun{
#' ## get the top diff. expressed genes, then plot them
#'
#' gcoreVlnPlot(myeloid_obj,
#'  genes = getTopNGenesSC(seurat_results$global_overall, min_logFC = 0.4, N = 6),
#'  grouping_var = 'stage',
#'  groups = c('Pre', 'Post'))
#' }
getTopNGenesSC <- function(result,
                           N = 30,
                           min_logFC = log2(1.5),
                           min_padj = 0.05,
                           min_pct = 0.25,
                           exclusion_patterns = c("^ENS", '^MT', '^RP[SL]'),
                           direction = "mixed",
                           lfc_col = 'avg_log2FC',
                           pval_col = 'p_val_adj')
{
  if (!is.data.frame(result)) {
    stop('Result is not a data.frame. Please provide output from Seurat::FindMarkers')
  }
  if (!is.null(exclusion_patterns)) {
    message('Excluding genes from return based on "exclusion_patterns" parameter. See help for details')
    for (pattern in exclusion_patterns) {
      result <- result[!grepl(pattern, rownames(result)),]
    }
  }

  filtered_results <- result %>%
    dplyr::filter(!is.na(.data$p_val_adj)) %>%
    dplyr::filter(.data$p_val_adj <= min_padj) %>%
    dplyr::filter(abs(.data$avg_log2FC) > min_logFC) %>%
    dplyr::filter(.data$`pct.1` > min_pct | .data$`pct.2` > min_pct) %>%
    tibble::rownames_to_column("Gene") %>%
    dplyr::arrange(.data$p_val_adj)

  if (direction == "mixed") {
    filtered_results <- filtered_results %>% utils::head(N)
  }
  else if (direction == "up") {
    filtered_results <- filtered_results %>%
      dplyr::filter(.data$avg_log2FC > 0) %>%
      utils::head(N)
  }
  else if (direction == "down") {
    filtered_results <- filtered_results %>%
      dplyr::filter(.data$avg_log2FC < 0) %>%
      utils::head(N)
  }
  else if (direction == "equal") {
    filtered_results_up <- filtered_results %>%
      dplyr::filter(.data$avg_log2FC < 0) %>% utils::head(round(N/2))
    filtered_results_down <- filtered_results %>%
      dplyr::filter(.data$avg_log2FC > 0) %>%
      utils::head(round(N/2))
    filtered_results <- rbind(filtered_results_up, filtered_results_down)
  }
  else {
    stop('Set "direction" to one of "equal", "mixed", "up", or "down"')
  }
  return(filtered_results$Gene)
}
