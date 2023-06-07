#' Generates diagnostic plots for evaluating Integration methods
#'
#' Generates (1) UMAP colored by sample label; (2) UMAP colored by cell labels; (3) UMAP colored by seurat_clusters (for all resolution specified); (4) A raster plot showing the number of cells of each given cell type assigned to each given seurat_clusters (for all resolution specified).
#'
#' @param plot_list list object to store plots in
#' @param seurat.obj Seurat object
#' @param subset_id Name of subset (e.g. "full", "myeloid", or "CD8-T")
#' @param integration_name Name of integration method (e.g. "merge", "STACAS", "Harmony")
#' @param sample_col Name of metadata column to use as sample name (e.g. "capID")
#' @param cell_labels Name of metadata column where cell annotation labels are stored
#' @param res Resolution(s) to use for `FindClusters()`. May be scalar or vector of multiple resolutions to loop over.
#'
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
#'        )
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
#' plot_list.annotation <- list()
#'
#' # Generate diagnostic plots
#' plot_list.annotation <-
#'   plotIntegrationDiagnostics(
#'     # global params
#'     plot_list = plot_list.annotation, seurat.obj = s.Harmony_v1,
#'     sample_col = "Sample_Name", subset_id = "full", cell_labels = "ImmGen.LCMV.collapsed.labels",
#'     res = seq(0.1, 1.0, 0.1),
#'     # this integration run
#'     integration_name = "Harmony_logNorm_split"
#'   )
#'
#' #' # Run a different version of the integration (harmony_norm_merged = TRUE)
#' s.Harmony_v2 <- runIntegration(s.split, integration_method = "harmony", norm_method = "logNorm",
#'                             group.by.vars = "Sample_Name", harmony_norm_merged = TRUE)
#'
#' # Generate diagnostic plots
#' plot_list.annotation <-
#'   plotIntegrationDiagnostics(
#'     # global params
#'     plot_list = plot_list.annotation, seurat.obj = s.Harmony_v2,
#'     sample_col = "Sample_Name", subset_id = "full", cell_labels = "ImmGen.LCMV.collapsed.labels",
#'     res = seq(0.1, 1.0, 0.1),
#'     # this integration run
#'     integration_name = "Harmony_logNorm_merged"
#'   )
#'
#' # Plot all the plots in 2x2 grids for Harmony_logNorm_split (Harmony_v1)
#' for (res in seq(0.1, 1.0, 0.1)) {
#'   p <- ggarrange(plot_list.annotation$full$Harmony_logNorm_split$Sample_Name,
#'                  plot_list.annotation$full$Harmony_logNorm_split$ImmGen.LCMV.collapsed.labels,
#'                  plot_list.annotation$full$Harmony_logNorm_split$seurat_clusters[[paste(res)]],
#'                  plot_list.annotation$full$Harmony_logNorm_split$clust.ann.compare$ImmGen.LCMV.collapsed.labels[[paste(res)]],
#'                  ncol=2, nrow = 2) %>%
#'     annotate_figure(top = paste0("res = ", res))
#'   print(p)
#' }
#'
#' #' # Plot all the plots in 2x2 grids for Harmony_logNorm_merged (Harmony_v2)
#' for (res in seq(0.1, 1.0, 0.1)) {
#'   p <- ggarrange(plot_list.annotation$full$Harmony_logNorm_merged$Sample_Name,
#'                  plot_list.annotation$full$Harmony_logNorm_merged$ImmGen.LCMV.collapsed.labels,
#'                  plot_list.annotation$full$Harmony_logNorm_merged$seurat_clusters[[paste(res)]],
#'                  plot_list.annotation$full$Harmony_logNorm_merged$clust.ann.compare$ImmGen.LCMV.collapsed.labels[[paste(res)]],
#'                  ncol=2, nrow = 2) %>%
#'     annotate_figure(top = paste0("res = ", res))
#'   print(p)
#' }
#'
#'}
#' @export
plotIntegrationDiagnostics <- function(plot_list, s, subset_id, integration_name, sample_col, cell_labels, res) {
  # Save plots for comparing integration methods
  print("Saving plots.")
  plot_list[[subset_id]][[integration_name]][[sample_col]] <-
    plot.umap.integrated(seurat.obj, group.by=sample_col,
                         title=paste0(integration_name,":\n", sample_col),
                         label.size=0, legend = TRUE)

  plot_list[[subset_id]][[integration_name]][[cell_labels]] <-
    plot.umap.integrated(seurat.obj, group.by=cell_labels,
                         title=paste0(integration_name,":\n", cell_labels),
                         label.size=2, legend = TRUE)

  # Save plots for comparing clustering to annotations
  for (resi in res) {
    print(resi)
    seurat.obj <- FindClusters(seurat.obj, resolution = resi, assay = seurat.obj@active.assay, verbose = F)

    plot_list[[subset_id]][[integration_name]][["seurat_clusters"]][[as.character(resi)]] <-
      plot.umap.integrated(seurat.obj, group.by="seurat_clusters",
                           title=paste0(integration_name,":\n", "seurat_clusters"),
                           label.size=3, legend = TRUE)

    plot_list[[subset_id]][[integration_name]][["clust.ann.compare"]][[cell_labels]][[as.character(resi)]] <-
      plot.cluster.annot.compare(seurat.obj, labels = cell_labels, res = resi)
  }

  plot_list[[subset_id]][[integration_name]][["clustree"]] <-
    clustree(seurat.obj, prefix = paste0(seurat.obj@active.assay,"_snn_res."), verbose = F)

  return(plot_list)
}

# Plot UMAP of integrated (or merged) data
plot.umap.integrated <- function(seurat.obj, group.by, title, label.size, legend=TRUE) {
  p <- DimPlot(seurat.obj, reduction = "umap",
               group.by = group.by,
               label = F) +
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

# Modify a UMAP to fit a ggarranged format better
plot_smaller <- function(p) {
  p +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          text = element_text(size = 2)) +
    theme(legend.key.size = unit(0.5, "line"),
          legend.text = element_text(size=4),
          legend.spacing.x = unit(0, "line")) +
    guides(color = guide_legend(
      ncol=1,
      override.aes=list(shape=15,
                        size=1)))
}

# Plot raster/tile heatmap shoing % cells in a cluster assigned to each cell annotation
plot.cluster.annot.compare <- function(obj.seurat, labels, res = 0.1, assay = "RNA") {
  obj.seurat <- FindClusters(obj.seurat, resolution = res, assay = assay)

  p <- obj.seurat@meta.data %>%
    group_by(seurat_clusters, .data[[labels]]) %>%
    summarize(n=n()) %>%
    ungroup() %>%
    group_by(.data[[labels]]) %>%
    mutate(colSum = sum(n)) %>%
    ungroup() %>%
    dplyr::filter(colSum > 10) %>%
    ggplot(data=., aes(x=seurat_clusters, y = .data[[labels]],
                       fill = log10(n+10))) +
    geom_tile() +
    scale_fill_viridis(discrete=FALSE) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size=8),
          axis.title = element_blank(),
          legend.position = "top",
          panel.background = element_rect(fill = "black"))
}
