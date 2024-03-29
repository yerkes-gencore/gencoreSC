% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotIntegrationDiagnostics.R
\name{plotIntegrationDiagnostics}
\alias{plotIntegrationDiagnostics}
\title{Generate diagnostic plots for evaluating Integration methods}
\usage{
plotIntegrationDiagnostics(
  plot_list,
  seurat.obj,
  subset_id,
  integration_name,
  sample_col,
  cell_labels,
  res,
  plot_proportions = TRUE,
  ...
)
}
\arguments{
\item{plot_list}{list object to store plots in}

\item{seurat.obj}{Seurat object}

\item{subset_id}{Name of subset (e.g. "full", "myeloid", or "CD8-T")}

\item{integration_name}{Name of integration method (e.g. "merge", "STACAS", "Harmony")}

\item{sample_col}{Name of metadata column to use as sample name (e.g. "capID")}

\item{cell_labels}{Name of metadata column where cell annotation labels are stored}

\item{res}{Resolution(s) to use for \code{\link[Seurat:FindClusters]{FindClusters()}}. May be scalar or vector of multiple resolutions to loop over.}

\item{plot_proportions}{Whether to plot the proprotion of cells assigned to each cell type in a given cluster for \code{plotClusterAnnotTile()}}

\item{\dots}{Additional arguments passed to \code{Seurat::DimPlot()}, useful for controlling rasterization or point size for large datasets}
}
\value{
A list of ggplot objects in a structure like so:

list(
subset_id = list(
integration_name = list(
sample_col,
cell_labels,
"seurat_clusters" = list(
res # a list of resolutions if res is a vector such as res=seq(0.1,1.0,1.0)
),
"clust.ann.compare" = list(
cell_labels = list(
res
)
),
"clustree"
)
)
)
}
\description{
Generates (1) UMAP colored by sample label; (2) UMAP colored by cell labels; (3) UMAP colored by seurat_clusters (for all resolution specified); (4) A raster plot showing the number of cells of each given cell type assigned to each given seurat_clusters (for all resolution specified); (5) A clustree looking across all resolutions specified.
}
\note{
Please do not try to save the output as an RDS. The list of list of list structure inflates the filesize dramatically and it will hang or create an impractically large RDS.
Instead, it is better practice to save all of the seurat objects you want to compare and use this function to re-plot the plots as needed for reports etc.
}
\examples{
\dontrun{

# Run integration (harmony_norm_merged = FALSE)
s.Harmony_v1 <- runIntegration(s.split, integration_method = "harmony", norm_method = "logNorm",
                            group.by.vars = "Sample_Name", harmony_norm_merged = FALSE)

# Intiialize an empty list to store plots in
plot_list <- list()

# Generate diagnostic plots
plot_list <-
  plotIntegrationDiagnostics(
    # global params
    plot_list = plot_list, seurat.obj = s.Harmony_v1,
    sample_col = "Sample_Name", subset_id = "full", cell_labels = "ImmGen.labels",
    res = seq(0.1, 1.0, 0.1),
    # this integration run
    integration_name = "Harmony_split"
  )

#' # Run a different version of the integration (harmony_norm_merged = TRUE)
s.Harmony_v2 <- runIntegration(s.split, integration_method = "harmony", norm_method = "logNorm",
                            group.by.vars = "Sample_Name", harmony_norm_merged = TRUE)

# Generate diagnostic plots
plot_list <-
  plotIntegrationDiagnostics(
    # global params
    plot_list = plot_list, seurat.obj = s.Harmony_v2,
    sample_col = "Sample_Name", subset_id = "full", cell_labels = "ImmGen.labels",
    res = seq(0.1, 1.0, 0.1),
    # this integration run
    integration_name = "Harmony_merged"
  )

# Plot all the plots in 2x2 grids for Harmony_split (Harmony_v1)
for (res in seq(0.1, 1.0, 0.1)) {
  p <- ggarrange(plot_list$full$Harmony_split$Sample_Name,
                 plot_list$full$Harmony_split$ImmGen.labels,
                 plot_list$full$Harmony_split$seurat_clusters[[paste(res)]],
                 plot_list$full$Harmony_split$clust.ann.compare$ImmGen.labels[[paste(res)]],
                 ncol=2, nrow = 2) \%>\%
    annotate_figure(top = paste0("res = ", res))
  print(p)
}

#' # Plot all the plots in 2x2 grids for Harmony_merged (Harmony_v2)
for (res in seq(0.1, 1.0, 0.1)) {
  p <- ggarrange(plot_list$full$Harmony_merged$Sample_Name,
                 plot_list$full$Harmony_merged$ImmGen.labels,
                 plot_list$full$Harmony_merged$seurat_clusters[[paste(res)]],
                 plot_list$full$Harmony_merged$clust.ann.compare$ImmGen.labels[[paste(res)]],
                 ncol=2, nrow = 2) \%>\%
    annotate_figure(top = paste0("res = ", res))
  print(p)
}

}
}
