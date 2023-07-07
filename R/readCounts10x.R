#' Read 10X count data
#'
#' Reads in 10X counts and metadata and returns Seurat obj
#'
#' @param filepath path to UMI count matrix (can be raw or unfiltered)
#' @param capID capture ID label, assigned to 'project' in seurat object metadata
#' @param min.cells Include features detected in at least this many cells. Will subset the counts matrix as well. To reintroduce excluded features, create a new object with a lower cutoff.
#' @param min.features Include cells where at least this many features are detected.
#' @param strip.suffix Whether to strip suffix from cellranger cell ID
#' @param format Input format. "tsv" (default) uses `Seurat::Read10X` to load in for the barcodes.tsv, features.tsv and matrix.mtx in the directory specified by the `filepat`h param. "h5" loads the `*feature_bc_matrix.h5` file specifed by the `filepath` param.
#'
#' @return An object of class Seurat
#'
#' @importFrom Seurat Read10X CreateSeuratObject CreateAssayObject
#'
#' @note For `format = "tsv"` (default), `filepath` must specify the path to the *directory* containing `barcodes.tsv`, `features.tsv`, and `matrix.mtx`; whereas for `format = h5`, `filepath` must specify the path to the `*feature_bc_matrix.h5` *file*.
#'
#' @export
readCounts10x <- function(capID,
                          filepath,
                          min.cells=0,
                          min.features=0,
                          strip.suffix=FALSE,
                          format = "tsv") {

  if (format == "tsv") {
    counts_in <- Seurat::Read10X(data.dir = filepath)
  } else if (format == "h5") {
    counts_in <- Seurat::Read10X_h5(filepath, use.names = TRUE, unique.features = TRUE)
  } else {
    errorCondition("Must specify whether input format is 'tsv' or 'h5'")
  }

  # Accommodate count matrices files with only RNA and those with multiple
  if ("Gene Expression" %in% names(counts_in)) {
    obj <- Seurat::CreateSeuratObject(counts_in$`Gene Expression`,
                                      project = capID,
                                      min.cells = min.cells,
                                      min.features = min.features,
                                      strip.suffix = strip.suffix)
    for (assay_name in names(counts_in)){
      if (assay_name != 'Gene Expression'){
        assay_tag <- ifelse(assay_name == "Antibody Capture", "ADT",
                            ifelse(assay_name == "Peaks", "ATAC",
                                   assay_name))

        obj[[assay_tag]] <- Seurat::CreateAssayObject(counts = counts_in[[assay_name]])
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
