#' Plot QC metadata across multiple captures or samples as density ridges plot
#'
#' Plots QC metadata for multiple captures or samples with ggplot2::geom_density_ridges
#'
#' @param metadata Seurat object metadata after running 'addQCmetrics()' and 'addQCfilter()'
#' @param cutoffs named list of cutoffs
#' @param split.by metadata column to split ridge plots by (i.e. usually capture ID or sample ID)
#'
#' @return an object of class ggarrange, which is a ggplot or a list of ggplot.
#'
#' @export
plotQC_ridges <- function(metadata, cutoffs=NULL, split.by) {

  # Visualize the number of cell counts per sample
  p.nCells <- metadata %>%
    ggplot(aes(y=.data[[split.by]], group=.data[[split.by]], fill=.data[[split.by]]), color=NA) +
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    theme(legend.position="none") +
    ggtitle("N Cells")

  p1 <- plotQC_ridges.helper(metadata, x="nUMI", split.by = split.by,
                             log10x=T, cutoffs=cutoffs)

  p2 <- plotQC_ridges.helper(metadata, x="nGene", split.by = split.by,
                             log10x=T, cutoffs=cutoffs)

  p3 <- plotQC_ridges.helper(metadata, x="log10GenesPerUMI", split.by = split.by,
                             cutoffs=cutoffs)

  p.mtRatio <- plotQC_ridges.helper(metadata, x="mitoRatio", split.by = split.by,
                                    cutoffs=cutoffs)

  p.rbRatio <- plotQC_ridges.helper(metadata, x="riboRatio", split.by = split.by,
                                    cutoffs=cutoffs)

  allCaps_QC_ridges <- ggpubr::ggarrange(p.nCells, p1, p2, p3, p.mtRatio, p.rbRatio, ncol=2, nrow=3)
  return(allCaps_QC_ridges)
}

plotQC_ridges.helper <- function(metadata, x, split.by, cutoffs=NULL, log10x=F) {

  p <- metadata %>%
    ggplot(aes(x=.data[[x]], y=.data[[split.by]], fill=.data[[split.by]]), color=NA) +
    ggridges::geom_density_ridges(size = 0, alpha = 0.8) +
    theme_classic() +
    theme(legend.position="none") +
    ggtitle(x) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    theme(legend.position="none",
          axis.text.y = element_blank(), axis.title.y = element_blank(),
          axis.ticks.y = element_blank(), axis.line.y = element_blank())

  # Add vert line for min cutoff (if provided)
  cutoff.min <- cutoffs[[paste0(x,".min")]]
  if (!is.null(cutoff.min)) {
    if (cutoff.min != 0) {
      p <- p +
        geom_vline(xintercept = cutoff.min)
    }
  }
  # Add vert line for max cutoff (if provided)
  cutoff.max <- cutoffs[[paste0(x,".max")]]
  if (!is.null(cutoff.max)) {
    if (cutoff.max != Inf) {
      p <- p +
        geom_vline(xintercept = cutoff.max)
    }
  }

  if(log10x) {
    p <- p + scale_x_log10()
  }
  return(p)
}
