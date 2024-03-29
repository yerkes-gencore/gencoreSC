% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{IdentName}
\alias{IdentName}
\title{Get Seurat object \code{Ident} name}
\usage{
IdentName(obj)
}
\arguments{
\item{obj}{Seurat object}
}
\value{
Name of metadata column used as the active.ident.
}
\description{
Gets the metadata column name of the current Ident of the Seurat object
}
\details{
Will error if the \code{Ident} isn't derived from a metadata column. Not sure why that would ever happen though...
}
\examples{
\dontrun{
# Get name of Ident
IdentName(obj)
}

}
