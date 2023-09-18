#' Violin plot of gene(s) expression
#'
#' Similar to `Seurat::VlnPlot` but more customizable. Pass a seurat object and
#'  genes to plot. Optionally subset the data by a meta.data variable defined
#'  with `subset_var`, selecting level(s) `subset`. Useful to plot expression
#'  in one subset of your data without having to make separate objects.
#'  Optionally group the data
#'  by meta.data variable `grouping_var`, selecting levels `groups`. Optionally
#'  filter zeros to focus on cells expressing gene.
#'
#'  Requires `ggforce` to be installed
#'
#' @param obj Seurat object
#' @param genes Genes to plot
#' @param subset_var  Optional. Column of `obj@meta.data` to subset on. Default
#'  is `seurat_clusters` so you could `subset` on a specific cluster.
#' @param subset Levels of `subset_var` to subset data to
#' @param grouping_var Optional. Column of `obj@meta.data` to group data by.
#' @param groups Levels of `grouping_var` to include in the plot. Also used
#'  to specify order of levels.
#' @param filter_zeros Remove 0s from plot: default `TRUE`.
#' @param assay Assay to pull expression data from, default `RNA`
#'
#' @returns A ggplot
#' @export
#'
#' @examples
#' \dontrun{
#' ## Plots only the 'Intermediate' cells (as labeled in 'coarse_labels' column)
#' ## Groups results by 'Pre' and 'Post', as labeled in the 'stage' column'
#'   gcoreVlnPlot(myeloid_obj,
#'     genes = c('ISG15', 'ISG20', 'CD14', 'LDHA', 'IFITM1'),
#'     subset = 'Intermediate', subset_var = 'coarse_labels',
#'     grouping_var = 'stage', groups = c('Pre', 'Post'))
#' }
gcoreVlnPlot <- function(obj,
                         genes,
                         assay = 'RNA',
                         subset = NULL,
                         subset_var = 'seurat_clusters',
                         grouping_var,
                         groups,
                         filter_zeros = TRUE){
  if (!is.null(subset)){
    obj <- obj[,obj@meta.data[[subset_var]] %in% subset]
  }
  mat_to_plot <- reshape2::melt(as.matrix(obj@assays[[assay]]@data[genes,]))
  mat_to_plot <- merge(mat_to_plot,
                       obj@meta.data %>%
                         as.data.frame() %>%
                         dplyr::select(.data[[grouping_var]]),
                       by.x = 'Var2', by.y = 0)
  mat_to_plot[[grouping_var]] <- factor(mat_to_plot[[grouping_var]],
                                        levels = groups)
  mat_to_plot <- mat_to_plot %>%
    filter(.data[[grouping_var]] %in% groups)

  if (filter_zeros) {mat_to_plot <- mat_to_plot %>% filter(.data[['value']] != 0)}

  ggplot2::ggplot(mat_to_plot,
                  aes(x = .data[[grouping_var]],
                      y = .data[['value']],
                      color = .data[[grouping_var]])) +
    # geom_jitter() +
    ggplot2::geom_violin(draw_quantiles = 0.5) +
    ggforce::geom_sina(size = 0.01, alpha = 0.2) +
    ggplot2::facet_wrap(~Var1, scale='free_y') +
    ggplot2::theme_bw() +
    ggplot2::labs(caption = paste0(if(!is.null(subset)) {paste0('Showing expression of ', subset, ' cells\n')},
                                   if(!is.null(filter_zeros)) {'Only showing cells with non-zero expression'}),
                  y = 'Expression')
}
