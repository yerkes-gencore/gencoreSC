#' Read 10X count data
#'
#' Reads in 10X counts and metadata and returns Seurat obj
#'
#' @param filepath path to UMI count matrix (can be raw or unfiltered)
#' @param min.cells minimum number of cells per feature for RNA assay
#' @param min.features minimum number of features per cell for RNA assay
#' @param capID capture ID label, assigned to 'project' in seurat object metadata
#'
#' @return Seurat object with the number of cells and features removed stored in 'nRemoved' in 'misc' slot
#'
#' @examples
#' s <- readBD(counts_csv = here("data/Combined_Exact-Cells_DBEC_MolsPerCell.csv"), metadata_csv = here("data/Exact-Cells_Sample_Tag_Calls.csv"))
#' @note
#' To do:
#'   - add option to not remove multiplets and undetermined droplets
#'   - currently only accepts RNA or RNA + ADT assay, add ability to handle arbitrary additional assays (e.g. ATAC)
#' @export
readCounts <- function(filepath, min.cells=100, min.features=10, capID) {

  counts_in <- Seurat::Read10X(data.dir = filepath)

  # Accommodate count matrices files with only RNA and those with both RNA and ADT
  # This is uglier than re-assigning the assay to a new mat obj but saves memory
  if ("Gene Expression" %in% names(counts_in)) {
    # Create Seurat obj with no filtering
    dims0 <- (Seurat::CreateSeuratObject(counts_in$`Gene Expression`,
                                         project = capID,
                                         min.cells=0, min.features=0))@assays$RNA@counts %>%
      dim()
    # Filter min.cells and min.features and report numbers removed
    obj <- Seurat::CreateSeuratObject(counts_in$`Gene Expression`,
                                      project = capID,
                                      min.cells=min.cells, min.features=min.features)

  } else if(is.null(names(counts_in))) {
    print("Only one assay in counts matrix file. Assuming it is Gene Expression.")
    # Create Seurat obj with no filtering
    dims0 <- (Seurat::CreateSeuratObject(counts_in,
                                         project = capID,
                                         min.cells=0, min.features=0))@assays$RNA@counts %>%
      dim()
    # Filter min.cells and min.features and report numbers removed
    obj <- Seurat::CreateSeuratObject(counts_in,
                                      project = capID,
                                      min.cells=min.cells, min.features=min.features)
    print(paste("RNA assay added to Seurat obj with ", length(rownames(obj[["RNA"]])), "features."))
  } else {
    errorCondition("No Gene Expression assay in counts matrix?")
  }

  ## Check if ADT data are included, add to Seurat obj if so
  if ("Antibody Capture" %in% names(counts_in)) {
    # Add ADT UMI matrix
    obj[["ADT"]] <- CreateAssayObject(counts = counts_in$`Antibody Capture`)
    print(paste("ADT assay added to Seurat obj with the following", length(rownames(obj[["ADT"]])), "features:"))
    print(rownames(obj[["ADT"]]))
  } else {
    warning("No Antibody Capture assay found. Proceeding without ADT data.")
  }

  show(obj)
  obj@misc[["readCounts_nGenesRemoved"]] <- (dims0 - dim(obj@assays$RNA@counts))[1]
  obj@misc[["readCounts_nCellsRemoved"]] <- (dims0 - dim(obj@assays$RNA@counts))[2]
  print(paste("Removing", obj@misc[["readCounts_nGenesRemoved"]], "genes and", obj@misc[["readCounts_nCellsRemoved"]], "cells."))

  return(obj)
}
