#' Read 10X count data
#'
#' Reads in 10X counts and metadata and returns Seurat obj
#'
#' @param filepath path to UMI count matrix (can be raw or unfiltered)
#' @param capID capture ID label, assigned to 'project' in seurat object metadata
#' @inheritParams SeuratObject::CreateSeuratObject
#'
#' @return An object of class Seurat
#'
#' @importFrom Seurat Read10X CreateSeuratObject CreateAssayObject
#'
#' @export
readCounts10x <- function(capID,
                          filepath,
                          min.cells=0,
                          min.features=0,
                          strip.suffix=TRUE) {

  counts_in <- Seurat::Read10X(data.dir = filepath)
  # Accommodate count matrices files with only RNA and those with multiple
  if ("Gene Expression" %in% names(counts_in)) {
    obj <- Seurat::CreateSeuratObject(counts_in$`Gene Expression`,
                                      project = capID,
                                      min.cells = min.cells,
                                      min.features = min.features,
                                      strip.suffix = strip.suffix)
    for (name in names(counts_in)){
      if (name != 'Gene Expression'){
        obj[[name]] <- Seurat::CreateAssayObject(counts = counts_in[[name]])
      }
    }

  } else if(is.null(names(counts_in))) {
    print("Only one assay in counts matrix file. Assuming it is Gene Expression.")
    obj <- Seurat::CreateSeuratObject(counts_in,
                                      project = capID,
                                      min.cells = min.cells,
                                      min.features = min.features,
                                      strip.suffix = strip.suffix)
  } else {
    errorCondition("No Gene Expression assay in counts matrix?")
  }
  show(obj)
  return(obj)
}
