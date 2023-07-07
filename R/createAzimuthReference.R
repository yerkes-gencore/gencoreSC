#' Create an Azimuth reference and save to file
#'
#' Create an Azimuth compatible reference from an existing Seurat object
#'  with a counts matrix and cell-type annotations in a metadata column. This
#'  wrapper will process the reference in a standard SCTransform workflow.
#'  The column of metadata with annotations should be specified to the `metadata`
#'  parameter. The reference will be saved to the specified output folder to
#'  be referenced in calls to `Azimuth::RunAzimuth`.
#'
#' @param ref A Seurat object with a counts matrix and cell-type annotations
#' @param metadata_column Columns of metadata to transfer onto query datasets
#' @param output_folder Where to save the Azimuth reference files for future calls
#'  of `RunAzimuth`
#' @param ndims Number of dimensions to use for PCA and UMAP
#' @param \dots
#'
#' @returns A Seurat object with AzimuthData stored in the tools slot for use with Azimuth
#' @export
#'
#' @note Code based on scripts from azimuth-references github, such as
#'  https://github.com/satijalab/azimuth-references/blob/master/human_motorcortex/scripts/export.R
#'
#' @import Seurat
#' @import dplyr
#'
#' @examples
#' \dontrun{
#'   ## Load reference data as Seurat object
#'   ref <- readRDS(here('saved_rds/tabula_sapiens_lymph_subset.Rds'))
#'   ## No processing needed, plug in
#'   azimuth_ref <- createAzimuthReference(ref,
#'                                         'cell_type',
#'                                         output_folder = here('reference/lymph_node_TS'))
#'
#'   ## Don't need the reference object for this, just the files
#'   combined.obj <- RunAzimuth(combined.obj,
#'                              reference = here('reference/lymph_node_TS'),
#'                              umap.name = 'refUMAP.lnTS',
#'                              assay = 'SCT')
#' }
createAzimuthReference <- function(ref,
                                   metadata_column,
                                   output_folder,
                                   ndims = 50,
                                   ...){
  if (!dir.exists(output_folder)){
    message('Creating output dir')
    dir.create(output_folder, recursive = TRUE)
  } else {
    files <- dir(output_folder)
    if ('idx.annoy' %in% files | 'ref.Rds' %in% files){
      stop('Reference files already exist at the destination provided by "output_folder".\n
            Please provide a different path')
    }
  }
  if (ndims < 50){
    stop('Azimuth requires references to be built with at least 50 dimensions')
  }

  ## in case the object is a subset, remove unwanted label factors
  ref@meta.data <- ref@meta.data %>% droplevels()

  ref <- Seurat::SCTransform(ref,
                             #method = "glmGamPoi",
                             new.assay.name = "SCT")
  ref <- Seurat::RunPCA(ref,
                        assay='SCT',
                        npcs = ndims)
  ## Return model for later map projections
  ref <- Seurat::RunUMAP(ref,
                         dims = 1:ndims,
                         reduction = 'pca',
                         return.model = TRUE,
                         umap.method = "uwot")
  ## Create azimuth compatible reference
  ## Can't specify namespace here due to weird stuff with a tools() call
  ## https://github.com/satijalab/azimuth/issues/155
  ## And can't do an import here cause I want to leave azimuth
  ## as a suggested package, so settling for 2 notes
  ref <- AzimuthReference(object = ref,
                          refUMAP = "umap",
                          refDR = "pca",
                          refAssay = "SCT",
                          dims = 1:ndims,
                          metadata = metadata_column,
                          reference.version = "1.0.0",
                          k.param = 31,
                          ...)

  ## Save azimuth compatible ref
  Seurat::SaveAnnoyIndex(object = ref[["refdr.annoy.neighbors"]],
                         file = file.path(output_folder,
                                          "idx.annoy"))
  saveRDS(object = ref,
          file = file.path(output_folder, "ref.Rds"))
  return(ref)
}
