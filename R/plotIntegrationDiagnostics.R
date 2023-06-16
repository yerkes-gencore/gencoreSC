#' Generate diagnostic plots for evaluating Integration methods
#'
#' Generates (1) UMAP colored by sample label; (2) UMAP colored by cell labels; (3) UMAP colored by seurat_clusters (for all resolution specified); (4) A raster plot showing the number of cells of each given cell type assigned to each given seurat_clusters (for all resolution specified); (5) A clustree looking across all resolutions specified.
#'
#' @param plot_list list object to store plots in
#' @param seurat.obj Seurat object
#' @param subset_id Name of subset (e.g. "full", "myeloid", or "CD8-T")
#' @param integration_name Name of integration method (e.g. "merge", "STACAS", "Harmony")
#' @param sample_col Name of metadata column to use as sample name (e.g. "capID")
#' @param cell_labels Name of metadata column where cell annotation labels are stored
#' @param res Resolution(s) to use for \code{\link[Seurat:FindClusters]{FindClusters()}}. May be scalar or vector of multiple resolutions to loop over.
#' @param \dots Additional arguments passed to `Seurat::DimPlot()`, useful for controlling rasterization or
#'  point size for large datasets
#' @returns A list of ggplot objects in a structure like so:
#'
#' list(
#'  subset_id = list(
#'    integration_name = list(
#'      sample_col,
#'      cell_labels,
#'      "seurat_clusters" = list(
#'        res # a list of resolutions if res is a vector such as res=seq(0.1,1.0,1.0)
#'        ),
#'      "clust.ann.compare" = list(
#'        cell_labels = list(
#'          res
#'          )
#'        ),
#'      "clustree"
#'      )
#'    )
#'  )
#'
#' @note Please do not try to save the output as an RDS. The list of list of list structure inflates the filesize dramatically and it will hang or create an impractically large RDS.
#' Instead, it is better practice to save all of the seurat objects you want to compare and use this function to re-plot the plots as needed for reports etc.
#'
#' @examples
#' \dontrun{
#'
#' # Run integration (harmony_norm_merged = FALSE)
#' s.Harmony_v1 <- runIntegration(s.split, integration_method = "harmony", norm_method = "logNorm",
#'                             group.by.vars = "Sample_Name", harmony_norm_merged = FALSE)
#'
#' # Intiialize an empty list to store plots in
#' plot_list <- list()
#'
#' # Generate diagnostic plots
#' plot_list <-
#'   plotIntegrationDiagnostics(
#'     # global params
#'     plot_list = plot_list, seurat.obj = s.Harmony_v1,
#'     sample_col = "Sample_Name", subset_id = "full", cell_labels = "ImmGen.labels",
#'     res = seq(0.1, 1.0, 0.1),
#'     # this integration run
#'     integration_name = "Harmony_split"
#'   )
#'
#' #' # Run a different version of the integration (harmony_norm_merged = TRUE)
#' s.Harmony_v2 <- runIntegration(s.split, integration_method = "harmony", norm_method = "logNorm",
#'                             group.by.vars = "Sample_Name", harmony_norm_merged = TRUE)
#'
#' # Generate diagnostic plots
#' plot_list <-
#'   plotIntegrationDiagnostics(
#'     # global params
#'     plot_list = plot_list, seurat.obj = s.Harmony_v2,
#'     sample_col = "Sample_Name", subset_id = "full", cell_labels = "ImmGen.labels",
#'     res = seq(0.1, 1.0, 0.1),
#'     # this integration run
#'     integration_name = "Harmony_merged"
#'   )
#'
#' # Plot all the plots in 2x2 grids for Harmony_split (Harmony_v1)
#' for (res in seq(0.1, 1.0, 0.1)) {
#'   p <- ggarrange(plot_list$full$Harmony_split$Sample_Name,
#'                  plot_list$full$Harmony_split$ImmGen.labels,
#'                  plot_list$full$Harmony_split$seurat_clusters[[paste(res)]],
#'                  plot_list$full$Harmony_split$clust.ann.compare$ImmGen.labels[[paste(res)]],
#'                  ncol=2, nrow = 2) %>%
#'     annotate_figure(top = paste0("res = ", res))
#'   print(p)
#' }
#'
#' #' # Plot all the plots in 2x2 grids for Harmony_merged (Harmony_v2)
#' for (res in seq(0.1, 1.0, 0.1)) {
#'   p <- ggarrange(plot_list$full$Harmony_merged$Sample_Name,
#'                  plot_list$full$Harmony_merged$ImmGen.labels,
#'                  plot_list$full$Harmony_merged$seurat_clusters[[paste(res)]],
#'                  plot_list$full$Harmony_merged$clust.ann.compare$ImmGen.labels[[paste(res)]],
#'                  ncol=2, nrow = 2) %>%
#'     annotate_figure(top = paste0("res = ", res))
#'   print(p)
#' }
#'
#'}
#' @export
plotIntegrationDiagnostics <- function(plot_list,
                                       seurat.obj,
                                       subset_id,
                                       integration_name,
                                       sample_col,
                                       cell_labels,
                                       res,
                                       ...) {
  # Save plots for comparing integration methods
  print("Saving plots.")
  plot_list[[subset_id]][[integration_name]][[sample_col]] <-
    plotUmapIntegrated(seurat.obj, group.by=sample_col,
                       title=paste0(integration_name,":\n", sample_col),
                       label.size=0, legend = TRUE, ...)

  plot_list[[subset_id]][[integration_name]][[cell_labels]] <-
    plotUmapIntegrated(seurat.obj, group.by=cell_labels,
                       title=paste0(integration_name,":\n", cell_labels),
                       label.size=2, legend = TRUE, ...)

  # Save plots for comparing clustering to annotations
  for (resi in res) {
    print(resi)
    seurat.obj <- FindClusters(seurat.obj, resolution = resi, assay = seurat.obj@active.assay, verbose = F)

    plot_list[[subset_id]][[integration_name]][["seurat_clusters"]][[as.character(resi)]] <-
      plotUmapIntegrated(seurat.obj, group.by="seurat_clusters",
                         title=paste0(integration_name,":\n", "seurat_clusters"),
                         label.size=3, legend = TRUE, ...)

    plot_list[[subset_id]][[integration_name]][["clust.ann.compare"]][[cell_labels]][[as.character(resi)]] <-
      plotClusterAnnotTile(seurat.obj, labels = cell_labels, res = resi)
  }
  tryCatch({
    plot_list[[subset_id]][[integration_name]][["clustree"]] <-
      clustree::clustree(seurat.obj, prefix = paste0(seurat.obj@active.assay,"_snn_res."), verbose = F)
  }, error = function(e) {
    warning('\nError generating clustree, possibly due to not having clustering for more than one resolution in the object')
  }, warning = function(w){
    warning('\nError generating clustree, possibly due to not having clustering for more than one resolution in the object')
  })
  return(plot_list)
}

#' Plot UMAP of integrated (or merged) data
#'
#' Wrapper for DimPlot to plot aspect.ratio=1 plot with some other convenient formatting. Used internally by \code{\link[gencoreSC:plotIntegrationDiagnostics]{PlotIntegrationDiagnostics()}}.
#'
#' @param seurat.obj Seurat object
#' @param group.by Metadata column to group (color) UMAP by
#' @param title Title of plot
#' @param label.size Cluster and legend text label size
#' @param legend Whether to plot legend
#' @param \dots Additional arguments passed to `Seurat::DimPlot()`, useful for controlling rasterization or
#'  point size for large datasets
#'
#' @returns ggplot object UMAP
#'
#' @export
plotUmapIntegrated <- function(seurat.obj, group.by, title, label.size, legend=TRUE, ...) {
  p <- DimPlot(seurat.obj,
               reduction = "umap",
               group.by = group.by,
               label = F, ...) +
    ggtitle(title) +
    theme(aspect.ratio = 1,
          plot.title = element_text(size=8),
          axis.title = element_text(size=8),
          axis.text = element_text(size=8))
  if (!legend) {
    p <- p + theme(legend.position = "none",
                   legend.text = element_text(size = label.size))
  }
  p <- LabelClusters(p, id = group.by, size = label.size, repel = T, max.overlaps=Inf) +
    theme(legend.key.size = unit(0.5, "line"),
          legend.text = element_text(size=3),
          legend.spacing.x = unit(0, "line")) +
    guides(color = guide_legend(ncol=1, override.aes=list(shape=15, size=1)))
  return(p)
}

#' Modify a UMAP to fit a ggarranged format better
#'
#' Removes axes and shrinks legend size of \code{\link[Seurat:DimPlot]{DimPlot()}} or \code{\link[gencoreSC:plotUmapIntegrated]{plotUmapIntegrated()}}.
#'
#' @param p ggplot object
#' @param discrete whether to use guide_legend (discrete variables) or guide_colorbar (continuous variables)
#' @param smaller_legend whether to plot a smaller legend
#'
#' @returns ggplot object UMAP with no axes and smaller legend
#'
#' @export
# Modify a UMAP to fit a ggarranged format better
plot_smaller <- function(p, discrete=T, smaller_legend=T) {
  p <- p +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          text = element_text(size = 2))
  if (smaller_legend) {
    p <- p +
      theme(legend.key.size = unit(0.5, "line"),
            legend.text = element_text(size=4),
            legend.spacing.x = unit(0, "line"))
  }
  if (discrete) {
    p <- p +
      guides(color = guide_legend(ncol=1, override.aes=list(shape=15, size=8),
                                  label.theme = element_text(size = 10)))
  } else if (!discrete) {
    p <- p +
      guides(color = guide_colorbar(label.position = "right",
                                    barheight = unit(8, "line"),
                                    barwidth = unit(1, "line"),
                                    label.theme = element_text(size = 10)))
  }
  return(p)
}

#' Plot comparing unsupervised clusters to cell annotations
#'
#' Plot raster/tile heatmap showing the log10(n+10) cells in a cluster assigned to each cell type. Useful for manually assigning cell type labels to clusters. Used internally by \code{\link[gencoreSC:plotIntegrationDiagnostics]{PlotIntegrationDiagnostics()}}.
#'
#' @param obj.seurat Seurat object
#' @param labels Metadata column name for cell type annotation labels in metadata
#' @param res Clustering resolution, passed to \code{\link[Seurat:FindClusters]{Seurat::FindClusters()}}
#' @param assay Assay to use (default "RNA")
#'
#' @returns ggplot object UMAP
#'
#' @importFrom viridis scale_fill_viridis
#' @note There definitely better ways to visualize this. At present, the color is on an absolute scale which makes cell type assignment to smaller clusters appear less confident.
#'
#' Perhaps a heatmap scaled by column (cluster) would be better?
#'
#' @export
plotClusterAnnotTile <- function(obj.seurat, labels, res = 0.1, assay = "RNA") {

  obj.seurat <- FindClusters(obj.seurat, resolution = res, assay = assay)

  p <- obj.seurat@meta.data %>%
    group_by(.data[["seurat_clusters"]], .data[[labels]]) %>%
    summarize(n=n()) %>%
    ungroup() %>%
    group_by(.data[[labels]]) %>%
    mutate(colSum = sum(n)) %>%
    ungroup() %>%
    dplyr::filter(.data[["colSum"]] > 10) %>%
    ggplot(data=., aes(x=.data[["seurat_clusters"]], y = .data[[labels]],
                       fill = log10(n+10))) +
    geom_tile() +
    viridis::scale_fill_viridis(discrete=FALSE) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size=8),
          axis.title = element_blank(),
          legend.position = "top",
          panel.background = element_rect(fill = "black"))
}
