#' Remove cells classified as doublets
#'
#' @param obj An Seurat object with a metadata field for doublet calls
#' @param classifications_col Column of doublet calls
#' @param doublet_class String for identifying doulbets in `classifications_col`
#'
#' @returns An object of class Seurat
#' @export
#'
#' @examples
#' \dontrun{
#'   ## with scDblFinder
#'   objs <- lapply(objs, run_scDblFinder)
#'   objs <- lapply(objs, removeDFDoublets)
#'   ## With DoubletFinder
#'   objs <- lapply(objs, runDoubletFinder, PCs = 10, sct = FALSE, cores = 1)
#'   objs <- lapply(objs, removeDFDoublets,
#'    classifications_col = 'DF_classifications', doublet_class = 'Doublet')
#' }
removeDFDoublets <- function(obj,
                             classifications_col = 'DF_classifications',
                             doublet_class = ''
){
  doub_idx <- obj[[classifications_col]] == 'Doublet'
  print(paste0(obj@project.name, ': Removing ', as.character(sum(doub_idx)), ' doublet GEMS'))
  obj <- obj[,!doub_idx]
  obj@meta.data <- droplevels(obj@meta.data)
  return(obj)
}
