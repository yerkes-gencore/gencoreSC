#' Plot QC nUMI vs nGene
#'
#' Plots QC metadata for multiple captures or samples with ggplot2::geom_density_ridges
#'
#' @param metadata Seurat object metadata after running 'addQCmetrics()' and 'addQCfilter()'
#' @param cutoffs named list of cutoffs
#' @param split_by metadata column to split ridge plots by (i.e. usually capture ID or sample ID)
#' @param color_by metric with to use for scale_color_gradient (i.e. "mitoRatio" or "riboRatio" or any other metric with limits=c(0,1))
#'
#' @return an object of class ggarrange, which is a ggplot or a list of ggplot.
#'
#' @export
plotQC_joint <- function(metadata, split_by="capID", cutoffs = NULL,
                         color_by="mitoRato") {
  # Visualize the correlation between genes detected and number of UMIs and
  # determine whether strong presence of cells with low numbers of genes/UMIs
  p <- metadata %>%
    ggplot(aes(x=.data$nUMI, y=.data$nGene, color=.data[[color_by]])) +
      geom_point(alpha=0.2) +
      scale_color_gradient(low = "goldenrod1", high = "red", limits=c(0,1)) +
      stat_smooth(method=lm) +
      scale_x_log10() +
      scale_y_log10() +
      theme_classic() +
      facet_wrap(~.data[[split_by]])

  if (!is.null(cutoffs)) {
    p <- p +
      geom_vline(xintercept = cutoffs$nUMI.min) +
      geom_hline(yintercept = cutoffs$nGene.min)
  }
  return(p)
}
