#' Read 10X count data
#'
#' Reads in 10X counts and metadata and returns Seurat obj
#'
#' @param filepath path to UMI count matrix (can be raw or unfiltered)
#' @param capID capture ID label, assigned to 'project' in seurat object metadata
#' @param min.cells Include features detected in at least this many cells.
#'  Will subset the counts matrix as well.
#'  To reintroduce excluded features, create a new object with a lower cutoff.
#'
#' @return An object of class Seurat
#'
#' @importFrom Seurat Read10X CreateSeuratObject
#'
#' @export
readCounts10x <- function(filepath,
                          capID,
                          min.cells=0,
                          min.features=0) {

  counts_in <- Seurat::Read10X(data.dir = filepath)
  # Accommodate count matrices files with only RNA and those with multiple
  if ("Gene Expression" %in% names(counts_in)) {
    obj <- Seurat::CreateSeuratObject(counts_in$`Gene Expression`,
                                            project = capID,
                                            min.cells = min.cells,
                                            min.features=min.features)

  } else if(is.null(names(counts_in))) {
    print("Only one assay in counts matrix file. Assuming it is Gene Expression.")
    obj <- Seurat::CreateSeuratObject(counts_in,
                                            project = capID,
                                            min.cells = min.cells,
                                            min.features = min.features)
  } else {
    errorCondition("No Gene Expression assay in counts matrix?")
  }
  show(obj)
  return(obj)
}
