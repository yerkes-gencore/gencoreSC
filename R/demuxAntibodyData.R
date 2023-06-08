#' Wrapper for calling HTODemux on antibody capture data in a Seurat object
#'
#' Provide a list of sample labels and hashtag oligos as a named list,
#'  where the name is the antibody capture feature and the value is the
#'  desired label. The HTO data will be extracted into it's own assay
#'  (in case there is other antibody data in the original assay), have the
#'  feature names converted to the desired labels, and run `Seurat::HTODemux()`.
#'
#' @param obj A Seurat object with antibody capture / hashtag oligo capture data
#' @param labels A named character vector, where values are the desired labels
#'  and the names are the HTO features as output by Cellranger
#' @param assay Name of assay to pull data from
#' @inheritParams Seurat::NormalizeData
#' @param ... Additional arguments to pass to HTODemux
#'
#' @inheritDotParams Seurat::HTODemux -object
#'
#' @inherit Seurat::HTODemux return
#'
#' @import Seurat
#' @importFrom plyr mapvalues
#' @export
#'
#' @examples
#' \dontrun{
#' ##> hash_labels
#  ##  Hash1Human-C0251-TotalSeqC Hash2Human-C0252-TotalSeqC Hash3Human-C0253-TotalSeqC Hash4Human-C0254-TotalSeqC Hash5Human-C0255-TotalSeqC
#  ##                     "RAo18"                    "RCn19"                    "ROi18"                    "RSi19"                    "RTa19"
#'   objs_filt <- lapply(objs_filt, demuxAntibodyData, labels = hash_labels)
#' }
demuxAntibodyData <- function(obj,
                              labels,
                              assay = 'Antibody.Capture',
                              normalization.method = 'CLR',
                              ...){
  hto_data <- obj[[assay]]@counts[names(labels),]
  rownames(hto_data) <- plyr::mapvalues(rownames(hto_data),
                                        from=names(labels),
                                        to=labels)
  obj[['HTO']] <- Seurat::CreateAssayObject(counts = hto_data)
  obj <- Seurat::NormalizeData(obj,
                               assay='HTO',
                               normalization.method = normalization.method)
  Seurat::HTODemux(obj, assay='HTO', ...)
}
