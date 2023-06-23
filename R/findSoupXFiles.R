#' Find SoupX input files in Cellranger output folder
#'
#' Cellranger has different outputs based on what is ran, so this function
#'  should help locate the necessary files for SoupX from the labrynth of output
#'  folders. It has not been robustly tested with multiple Cellranger output
#'  types, so it might not be general enough for other runs, but hopefully
#'  there is at least some consistency in outputs.
#'
#' @param base_path Path to cellranger outputs
#' @param filtered_mat_pattern Pattern to find filtered outputs via `dir()`
#' @param unfiltered_mat_pattern Pattern to find unfiltered outputs via `dir()`
#' @param clusters_pattern Pattern to find cellranger clustering info
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
findSoupXFiles <- function(base_path,
                           filtered_mat_pattern = 'filtered_feature_bc_matrix$',
                           unfiltered_mat_pattern = 'raw_feature_bc_matrix$',
                           clusters_pattern = 'clustering'){
  filtered_mat_path <- dir(file.path(base_path, 'outs'),
                           recursive = TRUE, pattern = filtered_mat_pattern,
                           include.dirs = TRUE, full.names = TRUE)

  unfiltered_mat_path <- dir(file.path(base_path, 'outs'),
                             recursive = TRUE, pattern = unfiltered_mat_pattern,
                             include.dirs = TRUE, full.names = TRUE)

  clusters_path <- file.path(dir(file.path(base_path, 'outs'),
                                 recursive = TRUE, pattern = clusters_pattern,
                                 include.dirs = TRUE, full.names = TRUE),
                             'graphclust/clusters.csv')
  return(list(unfiltered_mat_path = unfiltered_mat_path,
              filtered_mat_path = filtered_mat_path,
              clusters_path = clusters_path))
}
