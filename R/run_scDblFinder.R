#' wrapper for running scDblFinder
#'
#' Converts Seurat object to SCE for scDblFinder and returns
#'  a seurat object with the scDblFinder outputs as added metadata.
#'
#' @param obj A minimally processed seurat object
#' @inheritParams scDblFinder::scDblFinder
#' @inheritDotParams scDblFinder::scDblFinder
#' @returns A seurat object with additional metadata columns
#' @export
#'
#' @import scDblFinder
#' @importFrom Seurat as.SingleCellExperiment
#' @examples
#' \dontrun{
#'   ## If you want to multithread, pass BPPARAM
#'   objs <- lapply(objs, run_scDblFinder, BPPARAM = BioCParallel::MulticoreParam(10, RNGseed = 45))
#' }
run_scDblFinder <- function(obj,
                            knownDoublets=NULL,
                            knownUse=NULL,
                            nfeatures = 2000,
                            ...){
  if (!is.null(knownDoublets)) {
    if (is.null(knownUse)) {
      errorCondition('Specify how you want known doublets to be used with "knownUse"')
    }
  }
  obj.sce <- Seurat::as.SingleCellExperiment(obj)
  obj.sce <- scDblFinder::scDblFinder(obj.sce,
                                      knownDoublets = knownDoublets,
                                      knownUse = knownUse,
                                      nfeatures = nfeatures,
                                      ...)
  obj$scDblFinder.score <- obj.sce$scDblFinder.score
  obj$scDblFinder.class <- obj.sce$scDblFinder.class
  obj
}
