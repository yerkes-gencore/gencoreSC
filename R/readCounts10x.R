#' Read 10X count data
#'
#' Reads in 10X counts and metadata and returns Seurat obj
#'
#' @param filepath path to UMI count matrix (can be raw or unfiltered)
#' @param capID capture ID label, assigned to 'project' in seurat object metadata
#' @inheritParams SeuratObject
#'
#' @return An object of class Seurat
#'
#' @importFrom SeuratObject CreateSeuratObject
#' @importFrom Seurat Read10X
#'
#' @export
readCounts10x <- function(filepath,
                          capID,
                          min.cells=0,
                          min.features=0) {

  counts_in <- Seurat::Read10X(data.dir = filepath)
  # Accommodate count matrices files with only RNA and those with multiple
  if ("Gene Expression" %in% names(counts_in)) {
    obj <- SeuratObject::CreateSeuratObject(counts_in$`Gene Expression`,
                                            project = capID,
                                            min.cells = min.cells,
                                            min.features=min.features)

  } else if(is.null(names(counts_in))) {
    print("Only one assay in counts matrix file. Assuming it is Gene Expression.")
    obj <- SeuratObject::CreateSeuratObject(counts_in,
                                            project = capID,
                                            min.cells = min.cells,
                                            min.features = min.features)
  } else {
    errorCondition("No Gene Expression assay in counts matrix?")
  }
  show(obj)
  return(obj)
}
