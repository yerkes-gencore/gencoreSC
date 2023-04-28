#' Plot both ridges and joint qc plots as a ggarrange object
#'
#' Plot both ridges and joint qc plots as an object of class ggarrange
#'
#' @param obj.unfilt Unfiltered Seurat object after running 'addQCmetrics()' and 'addQCfilter()'
#' @param filtName metadata column to filter by for the plot
#' @param split_by metadata column to split ridge plots by (i.e. usually capture ID or sample ID)
#' @param color_by metric with to use for scale_color_gradient (i.e. "mitoRatio" or "riboRatio" or any other metric with limits=c(0,1))
#' @param cutoffs named list of cutoff values used in 'addQCmetrics()' for the given filtName
#' @param title main title
#'
#' @return an object of class ggarrange, which is a ggplot or a list of ggplot.
#'
#' @export
plotQC_arranged <- function(obj.unfilt, filtName, split_by, color_by="mitoRatio", cutoffs, title) {
  obj.filt <- obj.unfilt[,which(SeuratObject::FetchData(obj.unfilt, vars=filtName) == T)]

  p_ridges <- obj.filt %>%
    .@meta.data %>%
    plotQC_ridges(cutoffs=cutoffs, split_by=split_by)

  p_joint <- obj.filt %>%
    .@meta.data %>%
    plotQC_joint(cutoffs=cutoffs, split_by=split_by, color_by=color_by)

  p_arranged <- ggpubr::ggarrange(p_ridges, p_joint, ncol=2, nrow=1) %>%
    ggpubr::annotate_figure(top = title)
  return(p_arranged)
}
