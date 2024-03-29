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
#' @param facet_colors Whether to color facet labels by `split_by`. Requires package `ggh4x`.
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
                            y.text = FALSE,
                            facet_colors = F) {

  # if it's a split object, merge metadata and subset according to logical in filtName
  if (!is.list(obj)) {
    md <- obj@meta.data %>%
      dplyr::filter(.[[filtName]])
  } else if (is.list(obj)) {
    md <- lapply(obj, function(x) {x@meta.data}) %>%
      dplyr::bind_rows() %>%
      dplyr::filter(.[[filtName]])
  }

  if (mixed_sort == T) {
    md <- md %>%
      dplyr::mutate(!!split_by := factor(.[[split_by]], levels = gtools::mixedsort(unique(md[[split_by]]))))
    # sample_order <- unique(md[[split_by]]) %>% gtools::mixedsort()
    # md[[split_by]] <- factor(md[[split_by]], levels = sample_order)
  }

  p_ridges <- md %>%
    plotQC_ridges(cutoffs=cutoffs, split_by=split_by)

  p_joint <- md %>%
    plotQC_joint(cutoffs=cutoffs, split_by=split_by, color_by=color_by, facet_colors = facet_colors)

  p_arranged <- ggpubr::ggarrange(p_ridges, p_joint, ncol=2, nrow=1) %>%
    ggpubr::annotate_figure(top = ggpubr::text_grob(title, face = "bold", size = 14))
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
    dplyr::mutate(!!split_by := forcats::fct_rev(.[[split_by]]))

  # Visualize the number of cell counts per sample
  p.nCells <- metadata %>%
    ggplot(aes(y=.data[[split_by]], group=.data[[split_by]], fill=.data[[split_by]]), color=NA) +
    geom_bar() +
    theme_classic() +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    theme(legend.position="none") +
    ggtitle("N Cells")

  p1 <- plotQC_ridges.helper(metadata, x="nUMI", split_by = split_by,
                             assay = "RNA", log10x=T, cutoffs=cutoffs, y.text = y.text)

  p2 <- plotQC_ridges.helper(metadata, x="nGene", split_by = split_by,
                             assay = "RNA", log10x=T, cutoffs=cutoffs, y.text = y.text)

  p3 <- plotQC_ridges.helper(metadata, x="log10GenesPerUMI", split_by = split_by,
                             assay = "RNA", cutoffs=cutoffs, y.text = y.text)

  p.mtRatio <- plotQC_ridges.helper(metadata, x="mitoRatio", split_by = split_by,
                                    assay = "RNA", cutoffs=cutoffs, y.text = y.text)

  p.rbRatio <- plotQC_ridges.helper(metadata, x="riboRatio", split_by = split_by,
                                    assay = "RNA", cutoffs=cutoffs, y.text = y.text)

  allCaps_QC_ridges <- ggpubr::ggarrange(p.nCells, p1, p2, p3, p.mtRatio, p.rbRatio, ncol=2, nrow=3)
  return(allCaps_QC_ridges)
}

plotQC_ridges.helper <- function(metadata, x, split_by, assay = "ADT", cutoffs=NULL, log10x=F, y.text = FALSE, binwidth = 1) {

  # Base ggplot object
  p <- metadata %>%
    ggplot(aes(x=.data[[x]], y=.data[[split_by]], height = stat(density), fill=.data[[split_by]]), color=NA)

  # Use density if RNA, use histogram if ADTs
  if (assay == "RNA") {
    p <- p +
      ggridges::geom_density_ridges(stat = "density", size = 0, alpha = 0.8)
  } else if (assay == "ADT") {
    p <- p +
      ggridges::geom_density_ridges(stat = "binline", binwidth = 1, draw_baseline = FALSE, size = 0, alpha = 0.8)
  }
  p <- p +
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
      facet_wrap(~.data[[split_by]])
  } else {
    # define facet_wrap strip color as a ggplot aesthetic
    num_splits <- unique(metadata[[split_by]]) %>% length()
    strip_colors <- ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(
      fill = scales::hue_pal()(num_splits) %>% rev
      ))

    p <- p +
      ggh4x::facet_wrap2(~.data[[split_by]], strip = strip_colors)
  }


  if (!is.null(cutoffs)) {
    p <- p +
      geom_vline(xintercept = cutoffs$nUMI.min) +
      geom_hline(yintercept = cutoffs$nGene.min)
  }
  return(p)
}

#' Plot ADT QC metadata across multiple captures or samples as histogram
#'
#' Plots ADT QC metadata across multiple captures or samples as a histogram using `geom_density_ridges(stat = "binline")`
#'
#' @param metadata Seurat object metadata after running 'addQCmetrics()' and 'addQCfilter()'
#' @param cutoffs Named list of cutoffs
#' @param split_by Metadata column to split ridge plots by (i.e. usually capture ID or sample ID or hash ID)
#' @param y.text Whether to include labels for split_by on the y axis. Can set to `TRUE` when labels are very short.
#'
#' @return An object of class ggarrange, which is a ggplot or a list of ggplot.
#'
#' @export
plotQC_ADTbinline <- function(metadata, cutoffs=NULL, split_by, y.text = FALSE) {

  metadata <- metadata %>%
    dplyr::mutate(!!split_by := forcats::fct_rev(.[[split_by]]))

  # Visualize the number of cell counts per sample
  p.nCells <- metadata %>%
    ggplot(aes(y=.data[[split_by]], group=.data[[split_by]], fill=.data[[split_by]]), color=NA) +
    geom_bar() +
    theme_classic() +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    theme(legend.position="none") +
    ggtitle("N Cells")

  p1 <- plotQC_ridges.helper(metadata, x="nCount_ADT", split_by = split_by,
                             assay = "ADT", binwidth = 1,
                             log10x=F, cutoffs=NULL, y.text = y.text)

  p2 <- plotQC_ridges.helper(metadata, x="nFeature_ADT", split_by = split_by,
                             assay = "ADT", binwidth = 1,
                             log10x=F, cutoffs=NULL, y.text = y.text)

  p3 <- plotQC_ridges.helper(metadata, x="log10ADTPerUMI", split_by = split_by,
                             assay = "ADT", binwidth = 0.05,
                             log10x=F, cutoffs=NULL, y.text = y.text)

  allCaps_QC_ridges <- ggpubr::ggarrange(p.nCells, p1, p2, p3, ncol=3, nrow=1)
  return(allCaps_QC_ridges)
}
