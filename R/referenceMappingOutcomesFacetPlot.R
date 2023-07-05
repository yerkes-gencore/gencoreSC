#' Plot reference mapping calls and scores facetted by cluster
#'
#' Creates violin plots of cell annotation calls and confidence scores, facetted
#'  by query object clusters. Metadata is pulled from columns in the Seurat
#'  object as specified by the user, so it should be agnostic to SingleR, Azimuth,
#'  or other reference mapping functions that provide labels and scores. You
#'  can optionally specify a minimum proportion threshold to plot, which will
#'  remove data for cell-types appearing in less than X percent of a cluster.
#'  Order of annotation levels appearing on the X axis can be specified with
#'  `facet_order`.
#'
#' @param obj A seurat object with columns for cell-type annotations and annotation scores
#' @param label_column Column to pull annotations labels
#' @param label_score_column Column to pull annotation scores
#' @param clusters_column Column to pull cell clusters
#' @param min_proportion  Minimum proportion of a cluster an annotation must
#'  meet to be plotted
#' @param ncol  Number of columns for facetting
#' @param n_font_size Font size of cluster size labels
#' @param facet_order Order of annotation levels to show on X axis
#'
#' @returns A ggplot object
#'
#' @import ggplot2
#' @import dplyr
#'
#' @export
#'
#' @examples
#' \dontrun{
#' map_prediction_facetplot(
#' combined.obj,
#' label_column = 'predicted.cell_type.lymphTS',
#' label_score_column = 'predicted.cell_type.score.lymphTS',
#' min_freq = 0.05
#' )
#' }
referenceMappingOutcomesFacetPlot <- function(obj,
                                     label_column,
                                     label_score_column,
                                     clusters_column = 'seurat_clusters',
                                     #mapping_score_column,
                                     min_proportion = 0,
                                     ncol = 3,
                                     n_font_size = 2.5,
                                     facet_order = NULL){

  # Calculate proportion of calls in each cluster
  sample_sizes <- obj@meta.data %>%
    dplyr::select({{ clusters_column }}, {{ label_column }}, {{ label_score_column }}) %>%
    dplyr::group_by(dplyr::across(c({{ clusters_column }},{{ label_column }}))) %>%
    dplyr::summarise(n = n(), .groups = 'drop_last') %>%
    dplyr::mutate('freq' = .data[[ 'n' ]] / sum(.data[[ 'n' ]])) %>%
    dplyr::filter(.data[[ 'freq' ]] > min_proportion)

  data <- obj@meta.data %>%
    dplyr::select({{ clusters_column }},
      {{ label_column }},
      {{ label_score_column }})

  # Old, for median mapping scores, specific to Azimuth
  # med_map_scores <- data %>%
  #   dplyr::select(seurat_clusters, {{ mapping_score_column }}) %>%
  #   group_by(seurat_clusters) %>%
  #   summarise(med = median(.data[[mapping_score_column]])) %>%
  #   mutate(label = paste0('Cluster ', seurat_clusters, '\nMedian mapping score: ', round(med, 3)))
  # medmapscores <- as.character(med_map_scores$label)
  # names(medmapscores) <- as.character(med_map_scores$seurat_clusters)

  plot_data <- merge(sample_sizes, data,
                     by = c(clusters_column, label_column),
                     all.x = TRUE)

  ## Relevel annotation levels if provided
  plot_data[[clusters_column]] <-
    factor(plot_data[[clusters_column]],
           levels = if (is.null(facet_order)){
             sort(unique(plot_data[[clusters_column]]), decreasing = FALSE)
             } else {
               facet_order})

  ggplot2::ggplot(plot_data, aes(x = .data[[ label_column ]],
                        y = .data[[ label_score_column ]],
                        color = .data[[ label_column ]])) +
    ggplot2::geom_violin(draw_quantiles = c(0.5)) +
    ggplot2::geom_jitter(size=0.2, alpha=0.35) +
    ggplot2::facet_wrap(. ~ .data[[ clusters_column ]],
               ncol = ncol,
               # labeller = ggplot2::labeller(.data[[ clusters_column ]] = medmapscores)
               ) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggplot2::geom_text(aes(label=n, y = 0.15), angle = 30, size = n_font_size) +
    ggplot2::labs(x='Predicted cell type',
                  color = 'Predicted cell type',
                  y = 'Prediction score',
                  title='Reference based predictions by cluster',
                  caption = paste0(
                    "Lines in violin plots indicate the median.",
                    "\nNumbers below violins indicate the number of cells in that cluster with that label.",
                    # '\nMedian mapping score is a different metric than the prediction score mapped on the Y axis',
                    if (min_proportion > 0) {
                      paste0("\nCalls for less than ",
                             min_proportion*100,
                             "% of a cluster population are omitted for clarity.")
                      } else NULL ))
}
