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
#'
#' singler_out <- runSingleR(obj.seurat = pbmc.merged,
#'                           ref = pbmc3k.final@assays$RNA@data,
#'                           labels = pbmc3k.final$seurat_annotations,
#'                           de.method = 'wilcox',
#'                           meta.prefix = 'singleR.pbmc3k')
#' ## Spliting the list
#' pbmc.merged <- singler_out[[1]]
#' singler_out <- singler_out[[2]]
#'
#' ## Using the `%<-%` pipe saves typing and memory
#' library(zealot)
#' c(pbmc.merged, singler_out) %<-% runSingleR(obj.seurat = pbmc.merged,
#'                                             ref = pbmc3k.final@assays$RNA@data,
#'                                             labels = pbmc3k.final$seurat_annotations,
#'                                             de.method = 'wilcox',
#'                                             meta.prefix = 'singleR.pbmc3k')
#'
#' ## Using ImmGen data for mouse immune cells
#' # Prep / load reference
#' ref.ImmGen <- celldex::ImmGenData()
#'
#' Get cell ontology names and other info
#' cl <- ontoProc::getOnto('cellOnto')
#'
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
  singleR.out <- SingleR::SingleR(test = singleR.out$logcounts,
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

#' Convert cell ontology IDs to names in SingleR output
#'
#' Convert bare cell ontology IDs used in SingleR to human-readable names
#'
#' @param metadata data.frame with labels, pruned.labels and delta.next columns. Typically either:
#'  - Seurat object metadata (after running SingleR), or
#'  - SingleR output object
#' @param cl Cell ontology database
#'  - Can be btained via: `cl <- ontoProc::getOnto('cellOnto')`
#' @param old.prefix The prefix used for the original runSingleR call (e.g. "ImmGen.ont")
#' @param new.prefix The new prefix to use for the human readable column (e.g. "ImmGen.ont.name")
#'
#' @returns data.frame with three new columns: `{new.prefix}.labels`, `{new.prefix}.pruned.labels`, `{new.prefix}.delta.next`.
#'
#'@examples
#' \dontrun{
#' ## Run on bare ontology IDs, then convert to human readable cell type labels
#' # Get cell ontology names and other info
#' cl <- ontoProc::getOnto('cellOnto')
#'
#' # ImmGen.ont
#' c(obj.seurat, SingleR.ImmGen.ont) %<-% obj.seurat %>%
#'   Run.SingleR(obj.seurat=.,
#'               ref=ref.ImmGen,
#'               labels=ref.ImmGen$label.ont,
#'               de.method="classic",
#'               meta.prefix = "ImmGen.ont")
#'
#' # These labels aren't very informative:
#' SingleR::plotScoreHeatmap(SingleR.ImmGen.ont, show.pruned = TRUE, fontsize=4)
#'
#' # Convert bare cell ontology IDs to informative names in seurat object
#' obj.seurat@meta.data <- convertClID2Name(metadata = obj.seurat@meta.data, cl = cl,
#'                                      old.prefix = "ImmGen.ont", new.prefix = "ImmGen.ont.name")
#'
#' # Convert bare cell ontology IDs to informative names in SingleR output and visualize
#' SingleR.ImmGen.ont.names <- convertClID2Name(metadata = SingleR.ImmGen.ont, cl = cl,
#'                                          old.prefix = "ImmGen.ont", new.prefix = "ImmGen.ont.name")
#'
#' SingleR::plotScoreHeatmap(SingleR.ImmGen.ont.names, show.pruned = TRUE, fontsize=4)
#' }
#'
#' @export
# Run SingleR
# Convert bare cell ontology IDs to informative names
convertClID2Name <- function(metadata, cl, old.prefix = "ImmGen.ont", new.prefix = "ImmGen.ont.name") {
  op <- old.prefix
  np <- new.prefix
  md <- metadata
  md[[paste0(np, ".labels")]] <- unname(cl$name[unname(md[[paste0(op, ".labels")]])])
  md[[paste0(np, ".pruned.labels")]] <- unname(cl$name[unname(md[[paste0(op, ".pruned.labels")]])])
  md[[paste0(np, ".delta.next")]] <- md[[paste0(op, ".delta.next")]]
  metadata <- md
  return(metadata)
}
