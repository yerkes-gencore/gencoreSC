#' Wrapper for running singleR
#'
#' Maps labels from a reference dataset onto a query set using SingleR. The
#'  default uses `BioCparallel::MulticoreParam` to multi-thread the execution,
#'  but that may need to be disabled using the `num.threads` argument if your
#'  machine doesn't support that execution.
#'
#' @param obj.seurat A Seurat object. Can be a merged, non-integrated set
#'  of all samples, as SingleR can opperate on a per-cell level.
#' @inheritParams SingleR::SingleR
#' @param meta.prefix Character vector used to name outputs stored in the query object
#' @param de.method Arguments controlling the choice of marker genes used for annotation,
#'  see SingleR::trainSingleR()
#' @inheritDotParams SingleR::SingleR -BNPARAM -BPPARAM
#'
#' @returns A list containing the query object with reference annotations added as 3 metadata
#'  columns based on the `meta.prefix` and a data frame of the full SingleR results.
#' @export
#'
#' @import SingleR
#' @importFrom scater logNormCounts
#' @importFrom SummarizedExperiment assays
#' @importFrom Seurat as.SingleCellExperiment
#' @examples
#' \dontrun{
#' ## Install reference dataset
#' SeuratData::InstallData('pbmc3k')
#' singler_out <- runSingleR(obj.seurat = pbmc.merged,
#'                           ref = pbmc3k.final@assays$RNA@data,
#'                           labels = pbmc3k.final$seurat_annotations,
#'                           de.method = 'wilcox',
#'                           meta.prefix = 'singleR')
#' ## Spliting the list
#' write_rds(singler_out[[2]], 'saved_rds/singleR_out.Rds')
#' pbmc.merged <- singler_out[[1]]
#' singler_out <- singler_out[[2]]
#' write_rds(pbmc.merged, 'saved_rds/merged_post_singleR.Rds')
#' }
runSingleR <- function(obj.seurat,
                       ref,
                       labels,
                       de.method = "classic",
                       meta.prefix = 'singleR',
                       num.threads = BiocParallel::bpnworkers(BiocParallel::MulticoreParam()),
                       ...) {
  singleR.out <- obj.seurat %>%
    Seurat::as.SingleCellExperiment(assay = "RNA") %>%
    scater::logNormCounts() %>%
    SummarizedExperiment::assays() # %>%
    # .$logcounts %>%
    SingleR::SingleR(test = singleR.out$logcounts,
                     ref = ref,
                     labels = labels,
                     de.method = de.method,
                     num.threads = num.threads,
                     ...)

  obj.seurat@meta.data[[paste0(meta.prefix,".labels")]] <- singleR.out$labels
  obj.seurat@meta.data[[paste0(meta.prefix,".pruned.labels")]] <- singleR.out$pruned.labels
  obj.seurat@meta.data[[paste0(meta.prefix,".delta.next")]] <- singleR.out$delta.next

  return(list(obj.seurat, singleR.out))
}
