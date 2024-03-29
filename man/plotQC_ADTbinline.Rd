% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotQCRidgesJoint.R
\name{plotQC_ADTbinline}
\alias{plotQC_ADTbinline}
\title{Plot ADT QC metadata across multiple captures or samples as histogram}
\usage{
plotQC_ADTbinline(metadata, cutoffs = NULL, split_by, y.text = FALSE)
}
\arguments{
\item{metadata}{Seurat object metadata after running 'addQCmetrics()' and 'addQCfilter()'}

\item{cutoffs}{Named list of cutoffs}

\item{split_by}{Metadata column to split ridge plots by (i.e. usually capture ID or sample ID or hash ID)}

\item{y.text}{Whether to include labels for split_by on the y axis. Can set to \code{TRUE} when labels are very short.}
}
\value{
An object of class ggarrange, which is a ggplot or a list of ggplot.
}
\description{
Plots ADT QC metadata across multiple captures or samples as a histogram using \code{geom_density_ridges(stat = "binline")}
}
