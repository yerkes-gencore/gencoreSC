#' Extract the top N cell-type annotations by cluster
#'
#' Extract the top N cell-type annotations for each cluster by number of cells.
#'  For each annotation per cluster the median confidence score is reported.
#'  The output can optionally include the median mapping score from Azimuth outputs.
#'
#' @param obj A Seurat object with cell-type annotations and scores from something
#'  like Azimuth or SingleR
#' @param top_n How many annotations to report for each cluster
#' @param label_column Metadata column to extract cell-type labels from
#' @param label_score_column Metadata column to extract annotation scores from
#' @param mapping_score_column Optional, column to extract mapping scores from
#'  (from Azimuth)
#' @param clusters_column Column to extract clusters from
#'
#'
#' @import dplyr
#' @importFrom stats median
#'
#' @returns A dataframe of `top_n` rows per cluster
#' @export
#'
#' @examples
#' \dontrun{
#' extractRefMapTopResults(combined.obj,
#' label_column = 'predicted.cell_type.lymphTS',
#' label_score_column = 'predicted.cell_type.score.lymphTS',
#' mapping_score_column = 'mapping.score.lymphTS',
#' top_n = 2) %>%
#'   knitr::kable(digits = 3,
#'                col.names = c('Cluster',
#'                              'Predicted cell type',
#'                              'Median label confidence score',
#'                              'Median mapping confidence score',
#'                              'Proportion of cluster'),
#'                caption = 'Top 2 most common calls for each cluster') %>%
#'   kable_styling(full_width = F)
#'
#' }
extractTopRefMapResults <- function(obj,
                                    top_n,
                                    label_column,
                                    label_score_column,
                                    mapping_score_column = NULL,
                                    clusters_column = 'seurat_clusters'){
  obj@meta.data %>%
    dplyr::select({{ clusters_column }},
                  {{ label_column }},
                  {{ label_score_column }},
                  {{ mapping_score_column }}) %>%
    dplyr::group_by(pick({{ clusters_column }}, {{ label_column }})) %>%
    dplyr::summarise(median_label_score = stats::median(.data[[ label_score_column ]]),
                     median_mapping_score = if (!is.null(mapping_score_column)){
                       stats::median(.data [[ mapping_score_column ]])
                       } else NULL,
                     n = n(),
                     .groups = 'drop_last') %>%
    dplyr::mutate('freq' = .data[[ 'n' ]] / sum(.data [[ 'n' ]])) %>%
    dplyr::slice_max(.data[[ 'freq' ]], n = top_n) %>%
    dplyr::select(-.data[[ 'n' ]]) %>%
    dplyr::arrange(desc(.data[[ 'freq' ]]), .by_group = TRUE)
}
