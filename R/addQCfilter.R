#' Add QC filters
#'
#' Adds logical metadata column for a given set of arbitrary cutoffs and/or outliers
#'
#' @param obj Seurat object
#' @param split.by Metadata column to split by for outlier detection if obj is not already split
#' @param filterName Name of new metadata column with logical filter
#' @param by_threshold Logical to control whether to filter on thresholds defined in 'cutoffs'
#' @param cutoffs Named list of cutoffs
#' @param by_outlier Logical to control whether to filter outliers
#' @param nmads Median absolute deviation used to define outliers
#'
#'
#' @return Seurat object with new metadata column. If passed a split object, will return split object. Same for a merged/single capture object.
#'
#' @export
# calculate metadata filter column but don't filter yet
addQCfilter <- function(obj, split.by = NULL, filterName = "input",
                        by_threshold = FALSE, cutoffs = defaultCutoffs,
                        by_outlier = FALSE, nmads = 4) {

  if (is.list(obj)) {
    obj <- lapply(obj, FUN = function(x) {
      x <- addQCfilter.helper(x, filterName = filterName,
                              by_threshold = by_threshold, cutoffs = cutoffs,
                              by_outlier = by_outlier, nmads = nmads)
    })
  } else if (!is.list(obj)) {
    if (!by_outlier) {
      obj <- addQCfilter.helper(obj, filterName = filterName,
                                by_threshold = by_threshold, cutoffs = cutoffs,
                                by_outlier = by_outlier, nmads = nmads)
    } else if (is.null(split.by)) {
      # split the object
      obj <- Seurat::SplitObject(obj, split.by = split.by)
      # lapply over split object with addQCfilter.helper
      obj <- lapply(obj, FUN = function(x) {
        x <- addQCfilter.helper(x, filterName = filterName,
                                by_threshold = by_threshold, cutoffs = cutoffs,
                                by_outlier = by_outlier, nmads = nmads)
      })
      obj <- merge(x=obj[[1]], y=obj[2:length(obj)])
    } else if (!is.null(split.by)) {
      warning("Outlier detection on a merged object is not recommended. Either split object or specify split.by param.")
    }
  }
    return(obj)
  }

addQCfilter.helper <- function(obj, filterName = "input",
                               by_threshold = FALSE, cutoffs = defaultCutoffs,
                               by_outlier = FALSE, nmads = 4) {

    print(paste0("Creating new metadata column to store filter logical: ", filterName))

    # Create new metadata column with cell-level filtering logical
    if(by_threshold) {
      metadata <- obj@meta.data
      metadata <- metadata %>% dplyr::mutate(
        cutoffsOnly =
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
    } else {
      metadata <- obj@meta.data
      metadata <- metadata %>%
        dplyr::mutate(cutoffsOnly = TRUE)
      obj@meta.data <- metadata
    }
    # Create outlier metadata column if not already present
    if(! "outliers" %in% colnames(obj@meta.data)) {
      obj <- detect_outliers(obj, nmads)
    }
    # Create final filtering logical column based on whether filtering by outlier, cutoffs or both
    if (!by_outlier) {
      metadata <- obj@meta.data
      metadata <- metadata %>%
        dplyr::mutate(!!filterName := .data$cutoffsOnly) %>%
        dplyr::select(-.data$cutoffsOnly)
      obj@meta.data <- metadata
    } else {
      metadata <- obj@meta.data
      metadata <- metadata %>%
        # set logical to match only cells that pass the cutoffs AND are not outliers
        dplyr::mutate(!!filterName := .data$cutoffsOnly & !.data$outlier) #%>%
      obj@meta.data <- metadata
    }
  return(obj)
}

# Cell-level metrics
defaultCutoffs <- list(nUMI.max = Inf,
                       nUMI.min = 500,
                       nGene.max = Inf,
                       nGene.min = 250,
                       log10GenesPerUMI.max = Inf,
                       log10GenesPerUMI.min = 0.80,
                       mitoRatio.max = 0.2,
                       mitoRatio.min = 0,
                       riboRatio.max = 0.4,
                       riboRatio.min = 0)

### Outlier detection for a vector of values, based on median absolute deviation
is_outlier_mad <- function(data, nmads){
  data.mad <- stats::mad(data)
  outlier = (
    (data < stats::median(data) - nmads * data.mad) |
      (stats::median(data) + nmads * data.mad < data)
  )
  return(outlier)
}

detect_outliers <- function(obj, nmads=4){
  ## Log transforming data to normalize it
  outliers <- (is_outlier_mad(log1p(obj$nUMI), nmads) |
                   is_outlier_mad(log1p(obj$nGene), nmads) |
                   is_outlier_mad(log1p(obj$log10GenesPerUMI), nmads))
  obj$outlier <- outliers
  return(obj)
}
