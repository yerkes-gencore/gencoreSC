#' Plot split umap
#'
#' Plot a split umap from a split SeuratObject
#'
#' @param obj.split Unfiltered Seurat object after running 'addQCmetrics()' and 'addQCfilter()'
#' @param group.by metadata column to color cells by; passed to Seurat::DimPlot 'group.by'
#' @param title.prefix title prefix; name of obj (i.e. capture or sample ID) will be appended
#' @param title.size size of title; size of axis text will title.size-1
#' @param label.size size of DimPlot labels
#'
#' @return an object of class ggarrange, which is a ggplot or a list of ggplot.
#'
#' @export
plotUMAP_split <- function(obj.split, group.by="seurat_clusters", label.size=4, title.prefix="seurat_clusters", title.size=8) {
  umaps <- list()
  for(capID in names(obj.split)) {
    umaps[[capID]] <-
      Seurat::DimPlot(obj.split[[capID]], reduction = "umap",
                      group.by = group.by,
                      label = T, repel = T, label.size = label.size) +
      theme(legend.position = "none",
            aspect.ratio = 1,
            plot.title = element_text(size=title.size),
            axis.title = element_text(size=title.size-1),
            axis.text = element_text(size=title.size-1)) +
      ggtitle(paste(title.prefix, capID))
  }
  return(umaps)
}
