#' Remove genes present in low number of cells
#'
#' Remove rows from Seurat object for genes present in fewer than `min.cells`
#'  cells. This is normally done in `SeuratObject::CreateSeuratObject()`, but separating
#'  this out allows SoupX to be added to the workflow. After SoupX processing,
#'  this function can be ran to filter the object as you would normally.
#'
#' @param obj Seurat object
#' @param min.cells Minimum number of cells a gene must appear in to be kept
#' @param calculate_only If `TRUE` only show how many genes would be removed,
#'  but return the unfiltered object. If `FALSE`, return the filtered object
#'
#' @returns An object of class Seurat
#'
#' @importFrom Matrix rowSums
#' @export
#'
#' @examples
#' \dontrun{
#'   ## See the impact of filtering
#'   objs <- lapply(objs, filterLowGenes, min.cells = 100)
#'   ## apply filter
#'   objs <- lapply(objs, filterLowGenes, min.cells = 100, calculate_only = FALSE)
#' }
filterLowGenes <- function(obj,
                           min.cells,
                           calculate_only = TRUE){
  removed_genes <- (Matrix::rowSums(obj@assays$RNA@counts > 0) < min.cells)
  message(paste0(obj@project.name, ': ', sum(removed_genes),
                 ' genes found in fewer than ', min.cells, ' cells'))
  if (!calculate_only){
    message('Removing low count genes. New object:')
    show(obj[!removed_genes, ])
    return(obj[!removed_genes, ])
  } else {
    return(obj)
  }
}
