#' featureHeatmapByCluster
#'
#' Create a heatmap of feature expression by cluster. Specify 'assay' to
#'  pull either RNA or ADT data. If no features are provided, the first
#'  'n_features_to_plot' (default 50) will be used. This is to make it easy
#'  to quickly plot all ADT data, while also making this function usable
#'  with RNA data, but not breaking your computer by trying to accidentally
#'  render all RNA features.
#' @param obj A Seurat object
#' @param assay Assay to pull feature expression data from
#' @param features Optional: specify which features to plot. If none provided,
#'  defaults to plotting first `n_features_to_plot`
#' @param n_features_to_plot Number of features to plot if no features are provided.
#'  Default 50 to avoid accidentally attempting to plot all features.
#' @param cluster_column Metadata column containing cluster identities. Default
#'  is `seurat_clusters`
#' @param \dots Additional arguments passed to pheatmap::pheatmap
#'
#' @inheritParams pheatmap::pheatmap
#' @returns A ggplot converted from a pheatmap
#' @export
#'
#' @importFrom pheatmap pheatmap
#' @importFrom ggplotify as.ggplot
#'
#' @examples
#' \dontrun{
#'   featureHeatmapByCluster(obj.myeloids, assay = 'ADT', features = c('CD27.1', 'IgD', 'CD11c'))
#' }
featureHeatmapByCluster <- function(obj,
                                    assay,
                                    features = NULL,
                                    main = obj@project.name,
                                    cluster_rows = FALSE,
                                    cluster_cols = FALSE,
                                    n_features_to_plot = 50,
                                    cluster_column = 'seurat_clusters',
                                    ...) {
  if (!(assay %in% names(obj@assays))) {
    errorCondition('specified assay not found in provided object')
  }
  feats <- rownames(obj[[assay]])
  if (!is.null(features)) {
    feats <- features[features %in% feats]
    if (length(feats) < 1) {
      errorCondition('None of the requested feataures found in assay')
    }
  } else if (length(feats) > n_features_to_plot) {
    message(paste0('Only plotting ', n_features_to_plot, ' features, you can adjust this with the "n_features_to_plot" parameter'))
    feats <- feats[1:n_features_to_plot]
  }
  DefaultAssay(obj) <- assay
  obj <- obj[feats,]
  plot_data <- cbind(obj@meta.data, as.data.frame(t(as.matrix(obj[[assay]]@data)))) %>%
    dplyr::group_by(.data[[ cluster_column ]]) %>%
    dplyr::summarize_at(.vars = feats, .funs = median) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(cluster_column)

  heatmap <- pheatmap::pheatmap(t(plot_data),
                                color = viridis::viridis(25, option = "B"),
                                fonstize_row = 10,
                                border_color = NA,
                                legend = FALSE,
                                treeheight_row = 10,
                                treeheight_col = 10,
                                fontsize_col = 10,
                                main = main,
                                cluster_rows = cluster_rows,
                                cluster_cols = cluster_cols,
                                silent = TRUE,
                                ...) %>%
    ggplotify::as.ggplot()
  return(heatmap)
}
