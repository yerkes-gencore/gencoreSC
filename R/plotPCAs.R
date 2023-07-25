#' Plot PCAs for a list of categorical features
#'
#' Take a vector of features and plot PCA plots for each, returning a list.
#'
#' @param obj Seurat object
#' @param features Vector of features to to group.by and split.by in PCA plots
#'
#' @return List of plot objects
#' @export
plotPCAs <- function(obj, features = c("Phase", "mitoRatio", "riboRatio")) {
  p.list <- list()
  for (feature in features) {
    p.list[[feature]] <- Seurat::DimPlot(obj,
                                         reduction = "pca",
                                         group.by= feature,
                                         split.by = feature)
  }
  return(p.list)
}
