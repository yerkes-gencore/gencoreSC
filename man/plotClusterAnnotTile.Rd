% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotIntegrationDiagnostics.R
\name{plotClusterAnnotTile}
\alias{plotClusterAnnotTile}
\title{Plot comparing unsupervised clusters to cell annotations}
\usage{
plotClusterAnnotTile(
  obj.seurat,
  labels,
  assay = "RNA",
  plot_proportions = TRUE
)
}
\arguments{
\item{obj.seurat}{Seurat object}

\item{labels}{Metadata column name for cell type annotation labels in metadata}

\item{assay}{Assay to use (default "RNA")}

\item{plot_proportions}{Whether to plot proportions instead of raw counts, default TRUE}
}
\value{
ggplot object UMAP
}
\description{
Plot raster/tile heatmap showing the log10(n+10) cells in a cluster assigned to each cell type. Useful for manually assigning cell type labels to clusters. Used internally by \code{\link[gencoreSC:plotIntegrationDiagnostics]{PlotIntegrationDiagnostics()}}.
}
\note{
There definitely better ways to visualize this. At present, the color is on an absolute scale which makes cell type assignment to smaller clusters appear less confident.

Perhaps a heatmap scaled by column (cluster) would be better?
}
