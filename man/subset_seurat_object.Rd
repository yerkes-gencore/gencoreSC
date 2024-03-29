% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{subset_seurat_object}
\alias{subset_seurat_object}
\title{Subset a Seurat object on multiple variables and values}
\usage{
subset_seurat_object(obj, subset_vars, subsets)
}
\arguments{
\item{obj}{Seurat object}

\item{subset_vars}{Columns of object metadata to subset on}

\item{subsets}{Values of \code{subset_vars} to filter for. If multiple \code{subset_vars}
are provided, you can specify multiple \code{subsets} values for each var by
providing a list of character vectors. See examples.}
}
\value{
A filtered Seurat object
}
\description{
Provide one or more variables to subset on with one or more values to
filter for. This extends the base Seurat subsetting function to allow
multiple subsets to be supplied easily as function arguments. This is
intended for use within other functions, but will still work to return
a subsetted seurat object if you want to use it on it's own, although
Seurat's implementation may be preferred if you aren't constrained to
trying to do multiple subsets in a single function call.
}
\examples{
\dontrun{
## If you only have one value to subset for, the list convention isn't too picky
subset_seurat_object(obj = obj,
  subset_vars = 'cell_type',
  subsets = 'ISG-rich')

## If you want to subset for multiple values of a single subset_var, do this
subset_seurat_object(obj = obj,
  subset_vars = 'cell_type',
  subsets = list(c('ISG-rich', 'Other')))

## If you want to subset for multiple values of multiple subset_vars
subset_seurat_object(obj = obj,
  subset_vars = c('cell_type', 'timepoint'),
  subsets = list(c('ISG-rich', 'Other'), c('pre')))
}
}
