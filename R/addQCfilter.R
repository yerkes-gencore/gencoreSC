#' Add QC filters
#'
#' Adds logical metadata column for a given set of arbitrary cutoffs and/or outliers
#'
#' @param obj Seurat object
#' @param filterName Name of new metadata column with logical filter
#' @param cutoffs Named list of cutoffs with the following values:
#'  nUMI.max, nUMI.min, nGene.max, nGene.min, log10GenesPerUMI.max,
#'  log10GenesPerUMI.min, mitoRatio.max, mitoRatio.min, riboRatio.max, riboRatio.min
#'   You can provide a single list of values to be used for all captures (not recommended),
#'   or you can provide a list of lists, where each entry is thresholds for a single capture.
#'   See `generate_capture_QC_cutoffs()` for help making a cutoffs list.
#'
#'
#' @returns Seurat object(s) with a new metadata column indicating if the cell is within thresholds
#' @examples
#' \dontrun{
#'   ## where objs is a list of Seurat objects for each capture
#'   cutoffs <- lapply(objs, generate_cutoffs)
#'   objs <- mapply(addQCfilter, objs, 'outliers', cutoffs)
#' }
#'
#' @export

# calculate metadata filter column but don't filter yet
addQCfilter <- function(obj,
                        filterName = NULL,
                        cutoffs = NULL) {
  if (is.null(filterName)){
    stop('Use the "filterName" argument to specify a name/metadata column to store filtering outcomes')
  }
  if (is.null(cutoffs)){
    message('No cutoffs provided. Consider running `generate_capture_QC_cutoffs()` to tailor cutoffs to your captures, or
            manually provide cutoffs.\nUsing arbitrary default values to filter which may not be appropriate for your dataset.')
    cutoffs <- list(nUMI.max = Inf,
                    nUMI.min = 500,
                    nGene.max = Inf,
                    nGene.min = 250,
                    log10GenesPerUMI.max = Inf,
                    log10GenesPerUMI.min = 0.80,
                    mitoRatio.max = 0.2,
                    mitoRatio.min = 0,
                    riboRatio.max = 0.4,
                    riboRatio.min = 0)
  }
  ## For ridges plot, need a bool to filter cells on
  obj@meta.data$input <- TRUE

  ## filtering
  metadata <- obj@meta.data
  metadata <- metadata %>% dplyr::mutate(
    !!(filterName) :=
      .data$nUMI >= cutoffs$nUMI.min &
      .data$nUMI <= cutoffs$nUMI.max &
      .data$nGene >= cutoffs$nGene.min &
      .data$nGene <= cutoffs$nGene.max &
      .data$log10GenesPerUMI >= cutoffs$log10GenesPerUMI.min &
      .data$log10GenesPerUMI <= cutoffs$log10GenesPerUMI.max &
      .data$mitoRatio >= cutoffs$mitoRatio.min &
      .data$mitoRatio <= cutoffs$mitoRatio.max &
      .data$riboRatio >= cutoffs$riboRatio.min &
      .data$riboRatio <= cutoffs$riboRatio.max)
  obj@meta.data <- metadata
  return(obj)
}

#' generate_capture_QC_cutoffs
#'
#' Generates default cutoffs for all metrics of interest based on median
#' absolute deviation for the distribution of metric values in the capture.
#'
#' @param capture A seurat object
#' @param nmads Number of median absolute deviations, default 4
#' @param metrics Columns of the capture metadata to calculate thresholds for
#'
#' @returns a list of min and max thresholds for each of `metrics`
#' @export
#'
#' @examples
#' \dontrun{
#'   cutoffs <- generate_cutoffs(objs$Apr05_HH_5A)
#'   ## multiple captures
#'   cutoffs <- lapply(objs, generate_cutoffs)
#' }
generate_capture_QC_cutoffs <- function(capture,
                                        nmads = 4,
                                        metrics = c('nUMI',
                                                    'nGene',
                                                    'riboRatio',
                                                    'mitoRatio',
                                                    'log10GenesPerUMI')){
  cutoffs <- list()
  for (metric in metrics) {
    thresholds <- .generate_outlier_thresholds(capture, metric, nmads)
    cutoffs[[paste0(metric, '.max')]] <- thresholds[['upper']]
    cutoffs[[paste0(metric, '.min')]] <- thresholds[['lower']]
  }
  return(cutoffs)
}

## Generates upper and lower threshold based on MAD calculation
## helper function for generate_capture_QC_cutoffs
.generate_outlier_thresholds <- function(capture,
                                         metric,
                                         nmads,
                                         min.value = 0){
  data <- capture@meta.data[[metric]]
  data.mad <- stats::mad(data)
  data.med <- stats::median(data)
  upper <- data.med + nmads * data.mad
  lower <- data.med - nmads * data.mad
  if (lower < min.value) { lower <- min.value }
  return(list(upper = upper, lower = lower))
}
