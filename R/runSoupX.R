#' Wrapper for running SoupX
#'
#' Runs SoupX and saves the adjusted counts to a SoupChannel object for
#' downstream diagnostics and QC. Paths to input files can be provided manually,
#' or you can try to search for them with `findSoupXFiles()`.
#'
#' @param unfiltered_mat_path Path to unfiltered counts matrix file
#' @param filtered_mat_path   Path to Cellranger filtered counts matrix file
#' @param clusters_path       Path to Cellranger initial clustering data file
#' @inheritParams SoupX::autoEstCont
#'
#' @returns An object of class SoupChannel
#'
#' @import grDevices
#' @import SoupX
#' @importFrom Seurat Read10X
#' @export
#'
#' @examples
#' \dontrun{
#' soupx_files <- lapply(samples$FileID, function(x){
#'   findSoupXFiles(file.path(config$rootDir, config$alignmentDir, x))
#' })
#' names(soupx_files) <- samples$Label
#' sc <- lapply(soupx_files, function(x){
#'   runSoupX(x[['unfiltered_mat_path']], x[['filtered_mat_path']], x[['clusters_path']])
#' })
#' }
runSoupX <- function(unfiltered_mat_path,
                     filtered_mat_path,
                     clusters_path,
                     doPlot = TRUE){

  tod <- Seurat::Read10X(file.path(unfiltered_mat_path))

  if ("Gene Expression" %in% names(tod)) {
    tod <- tod$`Gene Expression`
  } else if(is.null(names(tod))) {
    tod <- tod
  } else {
    errorCondition("No Gene Expression assay in counts matrix?")
  }

  toc <- Seurat::Read10X(file.path(filtered_mat_path))
  if ("Gene Expression" %in% names(toc)) {
    toc <- toc$`Gene Expression`
  } else if(is.null(names(toc))) {
    toc <- toc
  } else {
    errorCondition("No Gene Expression assay in counts matrix?")
  }

  clus <- file.path(clusters_path)
  clus <- read.csv(clus)
  clusters <- clus$Cluster
  names(clusters) <- clus$Barcode

  sc <- SoupX::SoupChannel(tod = tod,
                    toc = toc)
  sc <- SoupX::setClusters(sc, clusters)
  sc <- SoupX::autoEstCont(sc, doPlot = doPlot)
  sc$adjusted_counts <- SoupX::adjustCounts(sc)
  sc$plot <- grDevices::recordPlot()
  sc
}
