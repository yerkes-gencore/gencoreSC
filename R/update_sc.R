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
