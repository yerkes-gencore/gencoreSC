% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotQueryRefMarkerHeatmaps.R
\name{plotQueryRefMarkerHeatmaps}
\alias{plotQueryRefMarkerHeatmaps}
\title{Plots heatmaps of the annotation marker gene expression}
\usage{
plotQueryRefMarkerHeatmaps(
  seurat_obj,
  ref_markers,
  labels,
  show_n_markers = 20,
  ...
)
}
\arguments{
\item{seurat_obj}{A seurat object (not a split object)}

\item{ref_markers}{Markers for the same labels in the reference dataset}

\item{labels}{Metadata column name for cell annotation labels}

\item{show_n_markers}{The number of markers to plot for each label. Default = 20.}

\item{...}{Arguments passed to scater::plotHeatmap}
}
\value{
A list of pheatmap objects
}
\description{
Plots a heatmap for each annotation label of the marker genes shared by the query and reference datasets. This is particularly useful for checking the markers determining label assignment by tools like SingleR.
}
\note{
Based on code from https://bioconductor.org/books/release/SingleRBook/annotation-diagnostics.html#based-on-marker-gene-expression
}
\examples{
\dontrun{
# Get markers for each annotation label
singlerMarkers_all <- metadata(singler_out)$de.genes

singlerMarker_heatmaps <- plotSingleRMarkerHeatmaps(seurat_obj = s,
                                                    ref_markers = singlerMarkers_all,
                                                    labels = "ts_heart_facs.labels",
                                                    show_n_markers = 30)
ggarrange(plotlist=singlerMarker_heatmaps, ncol = 3, nrow = 3)
}

}
