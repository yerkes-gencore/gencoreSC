#' Converts a Seurat object to SCE prepared for MAST
#'
#' Converts a processed Seurat object to a SingleCellExperiment object for
#'  use with MAST. Additional filtering is done on genes based on pattern expression
#'  to remove uninteresting genes, and filter for genes expressed in fewer than
#'  `min.cells` cells. This function also adds the cellular detection rate as
#'  metadata for modeling, per the MAST paper.
#'
#' @param obj A processed Seurat object with appropriate metadata for MAST modeling
#' @param exclusion_patterns Gene patterns to exclude from object, defaults to
#'  ribosomal and mitochondrial human patterns
#' @param min.cells Minimum number of cells a gene must appear in
#'
#' @returns An object of class SingleCellExperiment
#' @export
#'
#' @examples
#' \dontrun{
#' obj.sce <- prepare_MAST_obj(pdc_obj, min.cells = 150)
#' zlm_fit <- zlm(~ cdr + stage + Iga + stage:Iga + (1|individual),
#'   obj.sce,
#'   method = 'glmer',
#'   ebayes = FALSE,
#'   exprs_values = 'logcounts',
#'   parallel = TRUE,
#'   fitArgsD = list(nAGQ=0),
#'   silent = TRUE)
#' }
prepare_MAST_obj <- function(obj,
                             exclusion_patterns = c('^RP[SL]', '^MT'),
                             min.cells = 50){
  features <- rownames(obj)
  for (pattern in exclusion_patterns) {
    features <- features[!grepl(pattern, features)]
  }
  obj <- subset(obj, features = features)
  obj <- Seurat::as.SingleCellExperiment(obj, assay = 'RNA')

  ## Their thesholding model is based on binomial distribution, but I didn't observe that
  ## on our trimmed data. This may need to be explored more for generalization in the package.
  ## from vignette('MAITAnalysis', package = 'MAST')
  ## thresh <- MAST::thresholdSCRNACountMatrix(obj@assays@data$logcounts, nbins = 20, min_per_bin = 30)
  ## assays(sca, withDimnames = FALSE) <- list(thresh=thres$counts_threshold, tpm=assay(sca))

  ## For now, a simpler thresholding should suffice to weed out lowly expressed genes.

  features <- rowSums(obj@assays@data$counts > 0) > min.cells
  obj <- obj[features,]

  ## cellular detection rate
  obj$cdr <- colSums(SummarizedExperiment::assay(obj) > 0)
  obj$cdr <- scale(obj$cdr)
  obj
}


#' Run a Wald test on a MAST fit model
#'
#' Run a Wald test on a MAST fit model with a specified hypothesis and baseline
#'  to extract P value and LogFoldChange. The function is set up so the
#'  Wald test is run on the `hypothesis - baseline` and the logFC is calculated
#'  using the full `hypothesis` and `baseline` matrices.
#'
#'
#' @param zlmFit A MAST model fit produced from `MAST::zlm()`
#' @param hypothesis A character vector specifying the contrast of interest for the
#'  numerator. Must be compatible with `MAST::Hypothesis()`
#' @param baseline A MAST hypothesis object specifying the baseline/denominator
#'  for the comparison. Must be compatible with `MAST::Hypothesis()`
#' @param p_adj_method Method for multiple testing correction of pvalues
#'
#' @returns A dataframe of Wald test results
#' @export
#'
#' @examples
#' \dontrun{
#' ## These are based on the names of terms in my fit object
#' contrast1 <- '`(Intercept)`+`stagePost`'
#' contrast0 <- '`(Intercept)`'
#' res <- mastWaldTest(hypothesis = contrast1, baseline = contrast0, zlmFit = zlm_fit)
#' }
mastWaldTest <- function(zlmFit, hypothesis, baseline, p_adj_method = 'fdr'){
  result <- MAST::waldTest(zlmFit,
                           MAST::Hypothesis(paste0(hypothesis, '-(',baseline,')'), colnames(zlmFit@coefC)))
  result_lfc <- MAST::getLogFC(zlmFit,
                         contrast0 = MAST::Hypothesis(baseline, colnames(zlmFit@coefC)),
                         contrast1 = MAST::Hypothesis(hypothesis, colnames(zlmFit@coefC)))
  result <- merge(result[,,3], result_lfc, by.x = 0, by.y = 'primerid')
  result$padj <- p.adjust(result$hurdle, method = p_adj_method)
  result
}
