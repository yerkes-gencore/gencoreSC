#' Find SoupX input files in Cellranger output folder
#'
#' Cellranger has different outputs based on what is ran, so this function
#'  should help locate the necessary files for SoupX from the labrynth of output
#'  folders. It has not been robustly tested with multiple Cellranger output
#'  types, so it might not be general enough for other runs, but hopefully
#'  there is at least some consistency in outputs.
#'
#' @param base_path Path to cellranger outputs
#'
#' @returns A list of file paths
#' @export
#'
#' @examples
#' \dontrun{
#'
#' ## This returns all soupx files for all samples in the samplesheet
#' soupx_files <- lapply(samples$FileID, function(x){
#'   findSoupXFiles(file.path(config$rootDir, config$alignmentDir, x))
#' })
#' }
findSoupXFiles <- function(base_path){
  filtered_mat_path <- dir(file.path(base_path, 'outs'),
                           recursive = TRUE, pattern = 'filtered_feature_bc_matrix$',
                           include.dirs = TRUE, full.names = TRUE)

  unfiltered_mat_path <- dir(file.path(base_path, 'outs'),
                             recursive = TRUE, pattern = 'raw_feature_bc_matrix$',
                             include.dirs = TRUE, full.names = TRUE)

  clusters_path <- file.path(dir(file.path(base_path, 'outs'),
                                 recursive = TRUE, pattern = 'clustering', include.dirs = TRUE, full.names = TRUE),
                             'gene_expression_graphclust/clusters.csv')
  return(list(unfiltered_mat_path = unfiltered_mat_path,
              filtered_mat_path = filtered_mat_path,
              clusters_path = clusters_path))
}
