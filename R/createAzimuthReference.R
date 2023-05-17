#' Create an Azimuth reference and save to file
#'
#' Create an Azimuth compatible reference from an existing Seurat object
#'  with a counts matrix and cell-type annotations in a metadata column. This
#'  wrapper will process the reference in a standard SCTransform workflow.
#'  The column of metadata with annotations should be specified to the `metadata`
#'  parameter. The resulting
#'
#' @param ref A Seurat object with a counts matrix and cell-type annotations
#' @param metadata Columns of metadata to transfer onto query datasets
#' @param output_folder Where to save the Azimuth reference files for future calls
#'  of `Azimuth::RunAzimuth()`
#' @param dims Number of dimensions to use for PCA and UMAP
#' @param ... Additional arguments to pass to `Azimuth::AzimuthReference()`
#' @inheritDotParams Azimuth::AzimuthReference -refUMAP -refDR -object -refAssay
#'
#' @returns A Seurat object with AzimuthData stored in the tools slot for use with Azimuth
#' @export
#'
#' @import Azimuth
#' @import Seurat
#' @import dplyr
#'
#' @examples
#' \dontrun{
#'   ## Load reference data as Seurat object
#'   ref <- readRDS(here('saved_rds/tabula_sapiens_lymph_subset.Rds'))
#'   ## No processing needed, plug in
#'   azimuth_ref <- createAzimuthReference(ref, 'cell_type', output_folder = here('reference/lymph_node_TS'))
#'
#'   ## Don't need the reference object for this, just the files
#'   combined.obj <- RunAzimuth(combined.obj,
#'                              reference = here('reference/lymph_node_TS'),
#'                              umap.name = 'refUMAP.lnTS',
#'                              assay = 'SCT')
#' }
createAzimuthReference <- function(ref,
                                   metadata,
                                   output_folder,
                                   dims = 1:50,
                                   ...){

  files <- dir(output_folder)
  if ('idx.annoy' %in% files | 'ref.Rds' %in% files){
    stop('Reference files already exist at the destination provided by "output_folder".\n
                Please provide a different path')
  }

  ## in case the objectis a subset, remove unwanted label factors
  ref@meta.data <- ref@meta.data %>% droplevels()

  ref <- Seurat::SCTransform(ref)
  ref <- Seurat::RunPCA(ref)
  ## Return model for later map projections
  ref <- Seurat::RunUMAP(ref, dims = dims,
                 reduction = 'pca',
                 return.model = TRUE,
                 umap.method = "uwot")

  ## Create azimuth compatible reference
  ref <- Azimuth::AzimuthReference(object = ref,
                          refUMAP = "umap",
                          refDR = "pca",
                          refAssay = "SCT",
                          dims = dims,
                          metadata = metadata,
                          ...)

  ## Save azimuth compatible ref
  Azimuth::SaveAnnoyIndex(object = ref[["refdr.annoy.neighbors"]],
                 file = file.path(output_folder,
                                  "idx.annoy"))
  saveRDS(object = ref,
          file = file.path(output_folder, "ref.Rds"))
  return(ref)
}
