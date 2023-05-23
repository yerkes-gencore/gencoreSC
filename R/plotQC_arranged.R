#' Plot both ridges and joint qc plots as a ggarrange object
#'
#' Plot both ridges and joint qc plots as an object of class ggarrange
#'
#' @param obj Unfiltered Seurat object after running 'addQCmetrics()' and 'addQCfilter()'
#' @param filtName Metadata column to filter by for the plot
#' @param split_by Metadata column to split ridge plots by (i.e. usually capture ID or sample ID)
#' @param color_by Metric with to use for scale_color_gradient (i.e. "mitoRatio" or "riboRatio" or any other metric with limits=c(0,1))
#' @param cutoffs Named list of cutoff values used in 'addQCmetrics()' for the given filtName
#' @param title Main title
#' @param mixed_sort Whether to use `gtools::mixedsort()` to determine level order of `split_by`.
#' @param y.text Whether to include labels for split_by on the y axis on ridge plots. Can set to `TRUE` when labels are very short.
#'
#' @return An object of class ggarrange, which is a ggplot or a list of ggplot.
#'
#' @importFrom gtools mixedsort
#' @export
plotQCRidgesJoint <- function(obj,
                            filtName,
                            split_by,
                            color_by="mitoRatio",
                            cutoffs,
                            title,
                            mixed_sort = T,
                            y.text = FALSE) {

  # if it's a split object, merge metadata and subset according to logical in filtName
  if (!is.list(obj)) {
    md <- obj@meta.data %>%
      filter(.[[filtName]])
  } else if (is.list(obj)) {
    md <- lapply(obj, function(x) {x@meta.data}) %>%
      bind_rows() %>%
      filter(.[[filtName]])
  }

  if (mixed_sort == T) {
    md <- md %>%
      mutate(!!split_by := factor(.[[split_by]], levels = gtools::mixedsort(unique(md[[split_by]]))))
    # sample_order <- unique(md[[split_by]]) %>% gtools::mixedsort()
    # md[[split_by]] <- factor(md[[split_by]], levels = sample_order)
  }

  p_ridges <- md %>%
    plotQC_ridges(cutoffs=cutoffs, split_by=split_by)

  p_joint <- md %>%
    plotQC_joint(cutoffs=cutoffs, split_by=split_by, color_by=color_by)

  p_arranged <- ggpubr::ggarrange(p_ridges, p_joint, ncol=2, nrow=1) %>%
    ggpubr::annotate_figure(top = text_grob(title, face = "bold", size = 14))
  return(p_arranged)
}

#' Plot QC metadata across multiple captures or samples as density ridges plot
#'
#' Plots QC metadata for multiple captures or samples with ggplot2::geom_density_ridges
#'
#' @param metadata Seurat object metadata after running 'addQCmetrics()' and 'addQCfilter()'
#' @param cutoffs Named list of cutoffs
#' @param split_by Metadata column to split ridge plots by (i.e. usually capture ID or sample ID)
#' @param y.text Whether to include labels for split_by on the y axis. Can set to `TRUE` when labels are very short.
#'
#' @return An object of class ggarrange, which is a ggplot or a list of ggplot.
#'
#' @export
plotQC_ridges <- function(metadata, cutoffs=NULL, split_by, y.text = FALSE) {

  metadata <- metadata %>%
    mutate(!!split_by := forcats::fct_rev(.[[split_by]]))

  # Visualize the number of cell counts per sample
  p.nCells <- metadata %>%
    ggplot(aes(y=.data[[split_by]], group=.data[[split_by]], fill=.data[[split_by]]), color=NA) +
    geom_bar() +
    theme_classic() +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    theme(legend.position="none") +
    ggtitle("N Cells")

  p1 <- plotQC_ridges.helper(metadata, x="nUMI", split_by = split_by,
                             log10x=T, cutoffs=cutoffs, y.text = y.text)

  p2 <- plotQC_ridges.helper(metadata, x="nGene", split_by = split_by,
                             log10x=T, cutoffs=cutoffs, y.text = y.text)

  p3 <- plotQC_ridges.helper(metadata, x="log10GenesPerUMI", split_by = split_by,
                             cutoffs=cutoffs, y.text = y.text)

  p.mtRatio <- plotQC_ridges.helper(metadata, x="mitoRatio", split_by = split_by,
                                    cutoffs=cutoffs, y.text = y.text)

  p.rbRatio <- plotQC_ridges.helper(metadata, x="riboRatio", split_by = split_by,
                                    cutoffs=cutoffs, y.text = y.text)

  allCaps_QC_ridges <- ggpubr::ggarrange(p.nCells, p1, p2, p3, p.mtRatio, p.rbRatio, ncol=2, nrow=3)
  return(allCaps_QC_ridges)
}

plotQC_ridges.helper <- function(metadata, x, split_by, cutoffs=NULL, log10x=F, y.text = FALSE) {

  p <- metadata %>%
    ggplot(aes(x=.data[[x]], y=.data[[split_by]], fill=.data[[split_by]]), color=NA) +
    ggridges::geom_density_ridges(size = 0, alpha = 0.8) +
    theme_classic() +
    theme(legend.position="none") +
    ggtitle(x) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    theme(legend.position="none")

  if (y.text == FALSE) {
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
          axis.ticks.y = element_blank(), axis.line.y = element_blank())
  }

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

#' Plot QC nUMI vs nGene
#'
#' Visualize the correlation between genes detected and number of UMIs to determine whether strong presence of cells with low numbers of genes/UMIs
#'
#' @param metadata Seurat object metadata after running 'addQCmetrics()' and 'addQCfilter()'
#' @param cutoffs Named list of cutoffs
#' @param split_by Metadata column to split ridge plots by (i.e. usually capture ID or sample ID)
#' @param color_by Metric with to use for scale_color_gradient (i.e. "mitoRatio" or "riboRatio" or any other metric with limits=c(0,1))
#' @param facet_colors Logical; whether to plot facet colors with ggh4x::facet_wrap2 (experimental)
#'
#' @return An object of class ggarrange, which is a ggplot or a list of ggplot.
#'
#' @export
plotQC_joint <- function(metadata, split_by="capID", cutoffs = NULL,
                         color_by="mitoRato", facet_colors = FALSE) {

  metadata <- metadata %>%
    mutate(!!split_by := forcats::fct_rev(.[[split_by]]))

  # plot
  p <- metadata %>%
    ggplot(aes(x=.data$nUMI, y=.data$nGene, color=.data[[color_by]])) +
    geom_point(alpha=0.2) +
    scale_color_gradient(low = "goldenrod1", high = "red", limits=c(0,1)) +
    stat_smooth(method=lm) +
    scale_x_log10() +
    scale_y_log10() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  if (!facet_colors) {
    p <- p +
      facet_wrap(~fct_rev(.data[[split_by]]))
  } else {
    # define facet_wrap strip color as a ggplot aesthetic
    num_splits <- length(unique(metadata[[split_by]]))
    default_ggplot_hues <- scales::show_col(scales::hue_pal()(num_splits))
    strip_color <- ggh4x::strip_themed(
      background_x = ggh4x::elem_list_rect(fill = default_ggplot_hues)
    )

    p <- p +
      ggh4x::facet_wrap2(~fct_rev(.data[[split_by]]), strip = strip_color)
  }


  if (!is.null(cutoffs)) {
    p <- p +
      geom_vline(xintercept = cutoffs$nUMI.min) +
      geom_hline(yintercept = cutoffs$nGene.min)
  }
  return(p)
}

