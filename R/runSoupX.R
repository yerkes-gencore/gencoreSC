#' Wrapper for running SoupX
#'
#' Runs SoupX and saves the adjusted counts to a SoupChannel object for
#' downstream diagnostics and QC. Paths to input files can be provided manually,
#' or you can try to search for them with `findSoupXFiles()`.
#'
#' @param unfiltered_mat_path Path to unfiltered counts matrix file
#' @param filtered_mat_path   Path to Cellranger filtered counts matrix file
#' @param clusters_path       Path to Cellranger initial clustering data file
#' @param h5                  Whether the matrices are in h5 or mat format
#' @inheritParams SoupX::autoEstCont
#'
#' @returns An object of class SoupChannel
#'
#' @import grDevices
#' @import SoupX
#' @importFrom Seurat Read10X Read10X_h5
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
#'
#' ## Extract contamination estimates
#'
#' extractSoupXContamEst <- function(sc){
#'   rho <- sc$fit$rhoEst
#'   rho_low <- sc$fit$rhoFWHM[1]
#'   rho_high <- sc$fit$rhoFWHM[2]
#'   return(list(rho_low = rho_low, rho = rho, rho_high = rho_high))
#' }
#'
#' tmp <- data.table::rbindlist(lapply(sc, extractSoupXContamEst), idcol = '.id')
#' }
runSoupX <- function(unfiltered_mat_path,
                     filtered_mat_path,
                     clusters_path,
                     doPlot = TRUE,
                     h5 = TRUE){
  if (h5) {
    tod <- Seurat::Read10X_h5(file.path(unfiltered_mat_path))
  } else {
    tod <- Seurat::Read10X(file.path(unfiltered_mat_path))
  }


  if ("Gene Expression" %in% names(tod)) {
    tod <- tod$`Gene Expression`
  } else if(is.null(names(tod))) {
    tod <- tod
  } else {
    errorCondition("No Gene Expression assay in counts matrix?")
  }

  if (h5) {
    toc <- Seurat::Read10X_h5(file.path(filtered_mat_path))
  } else {
    toc <- Seurat::Read10X(file.path(filtered_mat_path))
  }
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

#' Find SoupX input files in Cellranger output folder
#'
#' Cellranger has different outputs based on what is ran, so this function
#'  should help locate the necessary files for SoupX from the labrynth of output
#'  folders. It has not been robustly tested with multiple Cellranger output
#'  types, so it might not be general enough for other runs, but hopefully
#'  there is at least some consistency in outputs.
#'
#' @param base_path Path to cellranger outputs
#' @param filtered_mat_pattern Pattern to find filtered outputs via `dir()`
#' @param unfiltered_mat_pattern Pattern to find unfiltered outputs via `dir()`
#' @param clusters_pattern Pattern to find cellranger clustering info
#'
#' @returns A list of file paths
#' @export
#'
#' @examples
#' \dontrun{
#'
#' ## This returns all soupx files for all samples in the samplesheet
#' soupx_files <- lapply(samples$FileID, function(x){
#'   findSoupXFiles(file.path(config$rootDir, config$alignmentDir, x))
#' })
#' }
findSoupXFiles <- function(base_path,
                           filtered_mat_pattern = 'filtered_feature_bc_matrix.h5$',
                           unfiltered_mat_pattern = 'raw_feature_bc_matrix.h5$',
                           clusters_pattern = 'clustering'){
  filtered_mat_path <- dir(file.path(base_path, 'outs'),
                           recursive = TRUE, pattern = filtered_mat_pattern,
                           include.dirs = TRUE, full.names = TRUE)

  unfiltered_mat_path <- dir(file.path(base_path, 'outs'),
                             recursive = TRUE, pattern = unfiltered_mat_pattern,
                             include.dirs = TRUE, full.names = TRUE)

  clusters_base_dir <- dir(file.path(base_path, 'outs'),
                           recursive = TRUE, pattern = clusters_pattern,
                           include.dirs = TRUE, full.names = TRUE)
  if (dir.exists(paste0(clusters_base_dir, "/graphclust"))) {
    clusters_path <- file.path(clusters_base_dir, "/graphclust/clusters.csv")
  } else if (dir.exists(paste0(clusters_base_dir, "/gene_expression_graphclust"))) {
    clusters_path <- file.path(clusters_base_dir, "/gene_expression_graphclust/clusters.csv")
  } else {
    errorCondition(paste0("Could not find cluster.csv file in ", clusters_base_dir, "/graphclust/ or gene_expression_graphclust/"))
  }

  return(list(unfiltered_mat_path = unfiltered_mat_path,
              filtered_mat_path = filtered_mat_path,
              clusters_path = clusters_path))
}

#' Updates a soup channel to only include cells remaining after processing and
#'  adds dimensional reduction loadings
#'
#' Takes a soup channel from a capture and a seurat object from that same
#'  capture to parse the soup channel based on what cells are left in the
#'  processed seurat object. This function also takes reduction loadings
#'  from the Seurat object and adds them to the soup channel, which
#'  allows additional plotting functions from `SoupX` to be used.
#'
#'  Since cell names may differ from the processed seurat object and soup channel,
#'  there is some functionality to make the names agree between the two
#'  objects. The defaults should catch most changes made by `Seurat::RenameCells()`
#'  and the suffixes add by `Seurat::CreateSeuratObject()`.
#'
#' @param sc A soup channel object
#' @param seur_obj A Seurat object processed to only contain cells of interest
#' @param reduction Reduction to use (default: 'umap')
#' @param gsub_start_pattern Pattern to replace at beginning of cell names
#' @param gsub_start_replace Replace term for `gsub_start_pattern`
#' @param gsub_end_pattern A second pattern to replace, after first is done.
#'  This is intended for trailing pattern, but technically is just another `gsub()` call.
#' @param gsub_end_replace Replace term for `gsub_end_pattern`
#' @param reduction Seurat reduction to pull cell embeddings from
#'
#' @importFrom SoupX setDR
#' @returns A soup channel object from `SoupX`
#' @export
#'
#' @examples
#' \dontrun{
#' ## Where 'sc' is a list of soup channels and 'objs_filt_mapped' is a list of seurat objects
#' sc_parsed <- mapply(update_sc, sc, objs_filt_mapped, SIMPLIFY = FALSE)
#'
#' ## soup channels can now be used with other SoupX plotting functions
#' plotMarkerMap(sc_parsed$R1, 'IGLV2-11')
#'
#' }
update_sc <- function(sc,
                      seur_obj,
                      reduction = 'umap',
                      gsub_start_pattern = '(^.+_)([ACTG]+.*)',
                      gsub_start_replace = '\\2',
                      gsub_end_pattern = '_\\d+$',
                      gsub_end_replace = ''){
  ## If cell names in your Seurat object are different than the
  ## original cell ranger outputs, you will need to adjust them
  ## to match. This uses gsub to trim a pattern from the start
  ## and end, assuming you added a prefix with RenameCells
  ## or had Seurat add the trailing number
  cells <- gsub(pattern = gsub_start_pattern,
                replacement = gsub_start_replace,
                x = colnames(seur_obj))
  cells <- gsub(pattern = gsub_end_pattern,
                replacement = gsub_end_replace,
                cells)
  sc$toc <- sc$toc[,cells]
  sc$metaData <- sc$metaData[cells,]
  sc$nDropUMIs <- sc$nDropUMIs[cells]
  sc$adjusted_counts <- sc$adjusted_counts[,cells]
  dr <- seur_obj@reductions[[reduction]]@cell.embeddings
  rownames(dr) <- cells
  sc <- SoupX::setDR(sc, dr)
  sc
}

#' Plot genes with high degree of contamination correction
#'
#' Uses the results of `runSoupX()` to produce a plot of genes with a high
#'  degree of contamination correction through SoupX. This can be defined as genes
#'  with a high portion of cell-level counts set to 0 after running SoupX (`plotSoupXGeneAdjustments()`).
#'  Another metric is the portion of reads removed for a given gene (`plotSoupXGeneAdjustments2()`).
#'  Ideally genes shown here should be strong cell-type markers,
#'  e.g. hemoglobins or immunoglobulins. For lightweight rendering
#'  and to focus on non-trivial signal, only genes found in a sufficient number
#'  of cells and with a sufficient degree of correction are plotted. You can
#'  also optionally just plot genes with a meaningful signal and exclude genes
#'  with only Ensembl IDs.
#'
#' @param sc SoupChannel object created from `runSoupX()`
#' @param ylim_min Minimum value of 'n'umber of original cells' to plot (y axis)
#' @param xlim_min Minimum value of 'portion of cells set to zero' to plot (x axis)
#' @param hide_ensembl Whether to include genes with `ens_pattern` in the plot
#' @param ens_pattern Regex pattern for identifying genes with no gene symbol.
#'  Only used to hide unknown genes from the plot if they aren't of interest
#'
#' @returns A plotly object
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr drop_na
#' @importFrom plotly ggplotly
#' @importFrom Matrix rowSums
#'
#' @examples
#' \dontrun{
#' sc <- runSoupX(tod = tod,
#'   toc = toc,
#'   clus = clusters)
#' soupXGeneAdjustPlot(sc, hide_ensembl = TRUE)
#' }
plotSoupXGeneAdjustments <- function(sc,
                                     ylim_min = 10,
                                     xlim_min = 0.05,
                                     hide_ensembl = FALSE,
                                     ens_pattern = '^ENS'){
  cntSoggy = Matrix::rowSums(sc$toc > 0)
  cntStrained = Matrix::rowSums(sc$adjusted_counts > 0)
  ratio = (cntSoggy - cntStrained) / cntSoggy
  count_diff <- data.frame(proportion = ratio,
                           original = Matrix::rowSums(sc$toc > 0)) %>%
    tibble::rownames_to_column('Gene') %>%
    tidyr::drop_na()

  if (hide_ensembl){
    count_diff <- count_diff[!grepl(ens_pattern, count_diff$Gene),]
  }
  soup_plot <- ggplot2::ggplot(count_diff,
                               aes(x=.data$proportion,
                                   y=.data$original,
                                   text=.data$Gene,
                                   alpha = 0.5)) +
    ggplot2::geom_point() +
    ggplot2::labs(x='Portion of cells set to zero',
                  y='Original number of cells with expression',
                  title = 'Genes most affected by SoupX correction') +
    ggplot2::theme_bw() +
    ggplot2::guides(alpha = 'none',
                    label = 'none') +
    ggplot2::lims(x = c(xlim_min, 1)) +
    ggplot2::scale_y_log10(limits=c((ylim_min),NA))

  message(paste0('This plot focuses on genes found in at least ',
                 ylim_min, ' cells in the original dataset',
                 'and a minimum adjustment of ',
                 xlim_min, ' to reduce rendering burden.'))
  plotly::ggplotly(soup_plot, tooltip = 'text') %>%
    plotly::config(
      displaylogo = FALSE,
      modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d", 'lasso2d', 'pan2d', 'autoScale2d', 'zoom2d')
    )
}


plotSoupXGeneAdjustments2 <- function(sc,
                                      ylim_min = 1e+02,
                                      xlim_min = 0.05,
                                      hide_ensembl = FALSE,
                                      ens_pattern = '^ENS'){
  compare <- 1 - (rowSums(sc$adjusted_counts) / rowSums(sc$toc))
  count_diff <- data.frame(diff=compare, total=rowSums(sc$toc)) %>%
    rownames_to_column('Gene') %>% drop_na()
  if (hide_ensembl){
    count_diff <- count_diff[!grepl(ens_pattern, count_diff$Gene),]
  }
  soup_plot <- ggplot2::ggplot(count_diff,
                               aes(x=.data$diff,
                                   y=.data$total,
                                   text=.data$Gene,
                                   alpha=0.5)) +
    ggplot2::geom_point() +
    ggplot2::labs(x='Portion of reads removed',
                  y='Original count total',
                  title = 'Genes most affected by SoupX correction',
                  caption = paste0('This plot focuses on genes with at least ',
                                   ylim_min, ' counts and a minimum adjustment of ',
                                   xlim_min, ' to reduce rendering burden.')) +
    ggplot2::theme_bw() +
    ggplot2::guides(alpha = 'none',
                    label = 'none') +
    ggplot2::lims(x = c(xlim_min, 1)) +
    ggplot2::scale_y_log10(limits=c((ylim_min),NA))
  plotly::ggplotly(soup_plot, tooltip = 'text') %>%
    plotly::config(
      displaylogo = FALSE,
      modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d", 'lasso2d', 'pan2d', 'autoScale2d', 'zoom2d')
    )
}

