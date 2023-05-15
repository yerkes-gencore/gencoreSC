#' Run SingleR
#'
#' Run SingleR::SingleR and return a seurat object with the new labels and scoring information in the metadata.
#'
#' @param obj.seurat Seurat object
#' @param ref A numeric matrix of (usually log-transformed) expression values from a reference dataset, or a SummarizedExperiment object containing such a matrix; see trainSingleR for details.
#' @param labels A character vector or factor of known labels for all samples in ref.
#' Alternatively, if ref is a list, labels should be a list of the same length. Each element should contain a character vector or factor specifying the label for the corresponding entry of ref.
#' @param de.method String specifying how DE genes should be detected between pairs of labels. Defaults to "classic", which sorts genes by the log-fold changes and takes the top de.n. Setting to "wilcox" or "t" will use Wilcoxon ranked sum test or Welch t-test between labels, respectively, and take the top de.n upregulated genes per comparison.
#' @param meta.prefix String appended to front of output metadata columns.
#'
#' @returns List with obj.seurat and SingleR output. Note: should the SingleR just be in misc. slot?
#'
#'@examples
#' \dontrun{
#' # Prep / load reference
#' ref.ImmGen <- celldex::ImmGenData()
#'
#' ## Run on highest level of reference
#' # ImmGen.main
#' list[obj.seurat, SingleR.ImmGen.main] <- obj.seurat %>%
#'   Run.SingleR(obj.seurat=.,
#'               ref=ref.ImmGen,
#'               labels=ref.ImmGen$label.main,
#'               de.method="classic",
#'               meta.prefix = "ImmGen.main")
#'
#' ## Run on bare ontology IDs, then convert to human readable cell type labels
#' # Get cell ontology names and other info
#' cl <- ontoProc::getOnto('cellOnto')
#'
#' # ImmGen.ont
#' list[obj.seurat, SingleR.ImmGen.ont] <- obj.seurat %>%
#'   Run.SingleR(obj.seurat=.,
#'               ref=ref.ImmGen,
#'               labels=ref.ImmGen$label.ont,
#'               de.method="classic",
#'               meta.prefix = "ImmGen.ont")
#'
#' # Convert bare cell ontology IDs to informative names
#' obj.seurat <- clID2clNames(obj.seurat, cl, old.prefix = "ImmGen.ont", new.prefix = "ImmGen.ont.name")
#' }
#'
#' @note
#' Use `de.method = "classic"` for bulk references (including microarray-based refs such as `celldex::ImmGenData()`).
#' Use `de.metohd = "wilcox"` for single cell references.
#'
#' @export
# Run SingleR
runSingleR <- function(obj.seurat, ref, labels, de.method, meta.prefix) {
  SingleR.out <- obj.seurat %>%
    as.SingleCellExperiment(assay= "RNA") %>%
    scater::logNormCounts() %>%
    assays() %>%
    .$logcounts %>%
    SingleR(test=.,
            ref=ref.ImmGen,
            labels=labels,
            de.method="classic")

  obj.seurat@meta.data[[paste0(meta.prefix,".labels")]] <- SingleR.out$labels
  obj.seurat@meta.data[[paste0(meta.prefix,".pruned.labels")]] <- SingleR.out$pruned.labels
  obj.seurat@meta.data[[paste0(meta.prefix,".delta.next")]] <- SingleR.out$delta.next

  return(list(obj.seurat, SingleR.out))
}

# Convert bare cell ontology IDs to informative names
clID2clNames <- function(obj.seurat, cl, old.prefix = "ImmGen.ont", new.prefix = "ImmGen.ont.name") {
  op <- old.prefix
  np <- new.prefix
  md <- obj.seurat@meta.data
  md[[paste0(np, ".labels")]] <- unname(cl$name[unname(md[[paste0(op, ".labels")]])])
  md[[paste0(np, ".pruned.labels")]] <- unname(cl$name[unname(md[[paste0(op, ".pruned.labels")]])])
  md[[paste0(np, ".delta.next")]] <- md[[paste0(op, ".delta.next")]]
  obj.seurat@meta.data <- md
  return(obj.seurat)
}
