#' Plot QC nUMI vs nGene
#'
#' Visualize the correlation between genes detected and number of UMIs to determine whether strong presence of cells with low numbers of genes/UMIs
#'
#' @param metadata Seurat object metadata after running 'addQCmetrics()' and 'addQCfilter()'
#' @param cutoffs named list of cutoffs
#' @param split_by metadata column to split ridge plots by (i.e. usually capture ID or sample ID)
#' @param color_by metric with to use for scale_color_gradient (i.e. "mitoRatio" or "riboRatio" or any other metric with limits=c(0,1))
#' @param facet_colors logical; whether to plot facet colors with ggh4x::facet_wrap2 (experimental)
#'
#' @return an object of class ggarrange, which is a ggplot or a list of ggplot.
#'
#' @export
plotQC_joint <- function(metadata, split_by="capID", cutoffs = NULL,
                         color_by="mitoRato", facet_colors = FALSE) {
  # plot
  p <- metadata %>%
    ggplot(aes(x=.data$nUMI, y=.data$nGene, color=.data[[color_by]])) +
      geom_point(alpha=0.2) +
      scale_color_gradient(low = "goldenrod1", high = "red", limits=c(0,1)) +
      stat_smooth(method=lm) +
      scale_x_log10() +
      scale_y_log10() +
      theme_classic()

  if (!facet_colors) {
    p <- p +
      facet_wrap(~.data[[split_by]])
  } else {
    # define facet_wrap strip color as a ggplot aesthetic
    num_splits <- length(unique(metadata[[split_by]]))
    default_ggplot_hues <- scales::show_col(scales::hue_pal()(num_splits))
    strip_color <- ggh4x::strip_themed(
      background_x = ggh4x::elem_list_rect(fill = default_ggplot_hues)
    )

    p <- p +
      ggh4x::facet_wrap2(~.data[[split_by]], strip = strip_color)
  }


  if (!is.null(cutoffs)) {
    p <- p +
      geom_vline(xintercept = cutoffs$nUMI.min) +
      geom_hline(yintercept = cutoffs$nGene.min)
  }
  return(p)
}
