#' Remove cells classified as doublets by DoubletFinder
#'
#' @param obj An Seurat object with a metadata field for DoubletFinder calls
#' @param classifications_col Column of DoubletFinder calls
#'
#' @returns An object of class Seurat
#' @export
#'
#' @examples
#' \dontrun{
#'   objs <- lapply(objs, runDoubletFinder, PCs = 10, sct = FALSE, cores = 1)
#'   objs <- lapply(objs, removeDFDoublets)
#' }
removeDFDoublets <- function(obj,
                             classifications_col = 'DF_classifications'
){
  doub_idx <- obj[[classifications_col]] == 'Doublet'
  print(paste0(obj@project.name, ': Removing ', as.character(sum(doub_idx)), ' doublet GEMS'))
  obj <- obj[,!doub_idx]
  obj@meta.data <- droplevels(obj@meta.data)
  return(obj)
}
