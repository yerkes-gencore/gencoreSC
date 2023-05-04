#' Small wrapper for running SoupX
#'
#' Runs SoupX and saves the adjusted counts to a SoupChannel object for
#' downstream diagnostics and QC.
#'
#' @inheritParams SoupX::SoupChannel
#' @inheritParams SoupX::setClusters
#'
#' @returns An object of class SoupChannel
#'
#' @import SoupX
#' @export
#'
#' @examples
#' \dontrun{
#' filt <- file.path(here::here(
#'   'cellranger_out/outs/per_sample_outs/sample/count/',
#'   'sample_filtered_feature_bc_matrix'))
#' filt <- Read10X(filt)
#' raw <- file.path(here::here(
#'   'cellranger_out/outs/per_sample_outs/sample/count/',
#'   'multi/count/raw_feature_bc_matrix'))
#' raw <- Read10X(raw)
#' clus <- file.path(here::here(
#'   '../cellranger_out/outs/per_sample_outs/sample/count/',
#'   'per_sample_outs/sample/count/',
#'   'analysis/clustering/gene_expression_graphclust/clusters.csv'))
#' clus <- read.csv(clus)
#' clusters <- clus$Cluster
#' names(clusters) <- clus$Barcode
#'
#' sc <- runSoupX(tod = raw$`Gene Expression`,
#'                toc = filt$`Gene Expression`,
#'                clus = clusters)
#' }
runSoupX <- function(tod,
                     toc,
                     clusters){
  sc <- SoupChannel(tod = tod,
                    toc = toc)
  sc <- setClusters(sc, clusters)
  sc <- autoEstCont(sc)
  sc$adjusted_counts <- adjustCounts(sc)
  sc
}
