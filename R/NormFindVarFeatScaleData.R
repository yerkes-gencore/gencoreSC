#' Normalize, FindVariableFeatures, ScaleData with optional blacklisting
#'
#' Normalize RNA counts with either SCTransform or log normalization; find variable features, with the option of blacklisting sets of genes; and opitonally scaling data.
#'
#' Each combination of options requires slightly different operations, so this wrapper helps keep code tidy and reproducible.
#'
#' @param s Seurat object
#' @param norm_method Normalization method
#' @param nfeatures Number of
regress.out = NULL,
feature.blacklist = NULL,
verbose = F,
scale_data = NULL
#'
#' @returns An object of class SoupChannel
#'
#' @import grDevices
#' @import SoupX
#' @importFrom Seurat Read10X
#' @export
#'
#' @examples
#' \dontrun{
#' soupx_files <- lapply(samples$FileID, function(x){
#'   findSoupXFiles(file.path(config$rootDir, config$alignmentDir, x))
#' })
#' names(soupx_files) <- samples$Label
#' sc <- lapply(soupx_files, function(x){
#'   runSoupX(x[['unfiltered_mat_path']], x[['filtered_mat_path']], x[['clusters_path']])
#' })
#' }
NormFindVarFeatScaleData <- function(s,
                                     norm_method = "logNorm",
                                     nfeatures = 2000,
                                     regress.out = NULL,
                                     feature.blacklist = NULL,
                                     verbose = F,
                                     scale_data = NULL) {
  # the native Seurat default for classic is to scale data, but for sctransform is to not scale data
  if (is.null(scale_data) & norm_method == "logNorm") {
    scale_data <- TRUE
  } else if (is.null(scale_data) & norm_method == "SCT") {
    scale_data <- FALSE
  }

  # scGate::genes.blacklist.default$Mm is a list of vectors, but we want a vector
  if (!is.null(feature.blacklist) & is.list(feature.blacklist)) {
    feature.blacklist <- unname(unlist(feature.blacklist))
  }
  # I don't think Seurat will try to normalize or scale SCT transformed data but just to be sure:
  DefaultAssay(s) <- "RNA"

  if (norm_method == 'logNorm'){
    s <- NormalizeData(s, verbose = verbose)
    s <- FindVariableFeatures(s,
                              selection.method = "vst",
                              nfeatures = nfeatures+length(feature.blacklist),
                              verbose = verbose)
    # Remove blacklisted features
    VariableFeatures(s) <- setdiff(VariableFeatures(s), feature.blacklist)
    # Remove extra features
    VariableFeatures(s) <- VariableFeatures(s)[1:nfeatures]
    if (scale_data == T) {
      print("scaling data")
      s <- ScaleData(s, verbose = verbose)
    }
  } else if (norm_method == 'SCT'){
    features <- setdiff(rownames(s), feature.blacklist)
    s <- SCTransform(s,
                     vst.flavor = "v2",
                     residual.features = features,
                     variable.features.n = nfeatures,
                     vars.to.regress = regress.out,
                     do.scale = scale_data,
                     verbose = verbose)
  } else if (!is.null(feature.blacklist) & !is.vector(feature.blacklist)) {
    errorCondition("feature.blacklist must be NULL, a vector or a list of vectors of gene symbols.")
  }
  return(s)
}
