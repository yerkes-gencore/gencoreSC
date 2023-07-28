#' Plots heatmaps of the annotation marker gene expression
#'
#' Plots a heatmap for each annotation label of the marker genes shared by the query and reference datasets. This is particularly useful for checking the markers determining label assignment by tools like SingleR.
#'
#' @param seurat_obj A seurat object (not a split object)
#' @param ref_markers Markers for the same labels in the reference dataset
#' @param labels Metadata column name for cell annotation labels
#' @param show_n_markers The number of markers to plot for each label. Default = 20.
#' @param ... Arguments passed to scater::plotHeatmap
#'
#' @returns A list of pheatmap objects
#'
#' @note
#' Based on code from https://bioconductor.org/books/release/SingleRBook/annotation-diagnostics.html#based-on-marker-gene-expression
#'
#' @examples
#' \dontrun{
#' # Get markers for each annotation label
#' singlerMarkers_all <- metadata(singler_out)$de.genes
#'
#' singlerMarker_heatmaps <- plotSingleRMarkerHeatmaps(seurat_obj = s,
#'                                                     ref_markers = singlerMarkers_all,
#'                                                     labels = "ts_heart_facs.labels",
#'                                                     show_n_markers = 30)
#' ggarrange(plotlist=singlerMarker_heatmaps, ncol = 3, nrow = 3)
#' }
#'
#' @export
plotQueryRefMarkerHeatmaps <- function(seurat_obj, ref_markers, labels, show_n_markers = 20, ...) {
  # Plotting aesthetics are much better with a shorter legend title
  seurat_obj@meta.data <- seurat_obj@meta.data %>%
    dplyr::mutate(labels = .data[[labels]])

  # Convert seurat obj to sce for input to scater::plotHeatmap
  obj_sce <- seurat_obj %>%
    Seurat::as.SingleCellExperiment(assay = "RNA") %>%
    scater::logNormCounts()

  # Find markers from the empirical data to compare to the reference based markers in `ref_markers`
  empirical.markers <- obj_sce %>%
    scran::findMarkers(., obj_sce[["labels"]], direction="up")

  # Plot a separate heatmap for each of the marker sets applied
  heatmaps <- lapply(names(ref_markers), function(x) {
    ref_markers_celltype <- unique(unlist(ref_markers[[x]]))
    m <- match(ref_markers_celltype, rownames(empirical.markers[[x]]))
    m <- ref_markers_celltype[rank(m) <= show_n_markers]

    heatmap <- scater::plotHeatmap(obj_sce, order_columns_by = "labels", features = m,
                                   silent = TRUE, main = x, ...)[[4]]
    return(heatmap)
  })
  names(heatmaps) <- names(ref_markers)
  return(heatmaps)
}
