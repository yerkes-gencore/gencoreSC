#' Normalize, FindVariableFeatures, ScaleData with optional blacklisting
#'
#' Normalize RNA counts with either SCTransform or log normalization; find variable features, with the option of blacklisting sets of genes; and opitonally scaling data.
#'
#' Each combination of options requires slightly different operations, so this wrapper helps keep code tidy and reproducible.
#'
#' @param s Seurat object
#' @param norm_method Normalization method (either "logNorm" = log transform or "SCT" = SCTransform)
#' @param nfeatures Number of features to pass to Seurat::FindVariableFeatures()
#' @param regress.out Features to regress out in SCTransform vars.to.regress
#' @param feature.blacklist Vector of genes to blacklist from variable features
#' @param verbose Print all messages from wrapped functions
#' @param scale_data Logical; whether to run Seurat::ScaleData(). If NULL (default), will scale for log normalization but not for SCTransform.
#'
#' @returns Seurat object
#' @export
#'
#' @examples
#' \dontrun{
#' # log normalize, find variable features, scale data
#' s <- NormFindVarFeatScaleData(s, norm_method = "logNorm")
#'
#' # SCTransform, find variable features, don't scale data
#' s <- NormFindVarFeatScaleData(s, norm_method = "SCT")
#'
#' # log normalize, find variable features blacklisting TCR gene (Mus musculus), scale data
#' library(scGate)
#' blacklist <- scGate::genes.blacklist.default$Mm$TCR
#' s <- NormFindVarFeatScaleData(s, norm_method = "logNorm", feature.blacklist = blacklist)
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
  Seurat::DefaultAssay(s) <- "RNA"

  if (norm_method == 'logNorm'){
    s <- Seurat::NormalizeData(s, verbose = verbose)
    s <- Seurat::FindVariableFeatures(s,
                                      selection.method = "vst",
                                      nfeatures = nfeatures+length(feature.blacklist),
                                      verbose = verbose)
    # Remove blacklisted features
    Seurat::VariableFeatures(s) <- setdiff(Seurat::VariableFeatures(s), feature.blacklist)
    # Remove extra features
    Seurat::VariableFeatures(s) <- Seurat::VariableFeatures(s)[1:nfeatures]
    if (scale_data == T) {
      print("scaling data")
      s <- Seurat::ScaleData(s, verbose = verbose)
    }
  } else if (norm_method == 'SCT'){
    features <- setdiff(rownames(s), feature.blacklist)
    s <- Seurat::SCTransform(s,
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
