#' Plot both ridges and joint qc plots as a ggarrange object
#'
#' Plot both ridges and joint qc plots as an object of class ggarrange
#'
#' @param obj Unfiltered Seurat object after running 'addQCmetrics()' and 'addQCfilter()'
#' @param filtName metadata column to filter by for the plot
#' @param split_by metadata column to split ridge plots by (i.e. usually capture ID or sample ID)
#' @param color_by metric with to use for scale_color_gradient (i.e. "mitoRatio" or "riboRatio" or any other metric with limits=c(0,1))
#' @param cutoffs named list of cutoff values used in 'addQCmetrics()' for the given filtName
#' @param title main title
#' @param sample_order Vector of levels for `split_by` used to order figures. Default is to use `gtools::mixedsort()` to deterimine levels.
#'
#' @return an object of class ggarrange, which is a ggplot or a list of ggplot.
#'
#' @importFrom gtools mixedsort
#' @export
plotQC_arranged <- function(obj,
                            filtName,
                            split_by,
                            color_by="mitoRatio",
                            cutoffs,
                            title,
                            sample_order = gtools::mixedsort(unique(obj[[split_by]]))[[1]]) {
  obj <- obj[,which(SeuratObject::FetchData(obj, vars=filtName) == T)]

  if (!is.null(sample_order)){
    if (all(sample_order %in% unique(obj[[split_by]][[1]]))){
      obj[[split_by]] <- factor(obj[[split_by]][[1]], levels = sample_order)
    } else {
      message("`sample_order` doesn't match levels of 'split_by' argument.\n
              Not able to order as requested.\n
              If you did not set sample_order, explicitly set it to NULL to avoid this error.")
    }
  }
  p_ridges <- obj %>%
    .@meta.data %>%
    plotQC_ridges(cutoffs=cutoffs, split_by=split_by)

  p_joint <- obj %>%
    .@meta.data %>%
    plotQC_joint(cutoffs=cutoffs, split_by=split_by, color_by=color_by)

  p_arranged <- ggpubr::ggarrange(p_ridges, p_joint, ncol=2, nrow=1) %>%
    ggpubr::annotate_figure(top = title)
  return(p_arranged)
}
