#' Plot split umap
#'
#' Plot a split umap from a split SeuratObject
#'
#' @param obj.split Unfiltered Seurat object after running \code{\link[gencoreSC:addQCmetrics]{addQCmetrics()}} and \code{\link[gencoreSC:addQCfilter]{addQCfilter()}}
#' @param reduction Reduction to use (default: "umap")
#' @param group.by Metadata column to color cells by; passed to \code{\link[Seurat:DimPlot]{DimPlot()}} 'group.by'
#' @param feature Metadata column to color cells by; passed to \code{\link[Seurat:FeaturePlot]{FeaturePlot()}} 'features' (but this function doesn't support a vector of features)
#' @param title.prefix Title prefix; name of obj (i.e. capture or sample ID) will be appended
#' @param label.size Size of DimPlot labels
#' @param plot_smaller Whether to plot smaller
#' @param legend Whether to plot legend (helpful to remove legend if plotting smaller)
#' @param ... Arguments passed to \code{\link[gencoreSC:plot_smaller]{plot_smaller()}}
#'
#' @return an object of class ggarrange, which is a ggplot or a list of ggplot.
#'
#' @export
# Plot UMAPs for split objects
plotUMAP_split <- function(obj.split, reduction = "umap", group.by=NULL, feature=NULL, title.prefix=NULL, label.size=4, plot_smaller=T, legend=F, ...) {
  umaps <- list()
  for(Sample_Name in names(obj.split)) {
    if (is.null(feature)) {
      umaps[[Sample_Name]] <- DimPlot(obj.split[[Sample_Name]], reduction = reduction,
                                      group.by = group.by,
                                      label = T, repel = T, label.size = label.size)
    } else if (!is.null(feature)) {
      umaps[[Sample_Name]] <- FeaturePlot(obj.split[[Sample_Name]], reduction = reduction,
                                          features = feature)#,
      #cols = viridis::viridis(25, option = "D"))
    }
    umaps[[Sample_Name]] <- umaps[[Sample_Name]] +
      theme(aspect.ratio = 1,
            plot.title = element_text(size=8),
            axis.title = element_text(size=8),
            axis.text = element_text(size=8)) +
      ggtitle(paste(title.prefix, Sample_Name))

    if (plot_smaller) {
      umaps[[Sample_Name]] <- plot_smaller(umaps[[Sample_Name]],
                                           discrete = is.null(feature),
                                           ...)
    }

    if (!legend) {
      umaps[[Sample_Name]] <- umaps[[Sample_Name]] +
        theme(legend.position = "none")
    }
  }
  return(umaps)
}
