#' A wrapper for DoubletFinder functions
#'
#' Takes a pre-processes Seurat object for a single capture and performs
#' DoubletFinder functions. The object is returned with an additional metadata
#' column of DoubletFinder classification calls. This can optionally be used with
#' Antibody capture data for hashed studies to provide 'ground truth' calls
#' and improve the doublet detection algorithm.
#'
#' @param obj A Seurat object for a single capture
#' @param PCs Number of PCs to use
#' @param sct Is data SCTransformed?
#' @param cores Number of cores to use for parallelization
#' @param ground_truth Optional: An nCell-length character vector of ground-truth doublet classifications (e.g., "Singlet" or "Doublet") used to gauge performance of logistic regression models trained using pANN vectors during ROC analysis.
#' @param doublet_rate Estimated doublet rate, based on sequencing device specifications, library prep, etc.
#' @param assay Assay to perform DoubletFinder on, default 'RNA'
#'
#' @returns A Seurat object with additional metadata columns, including doublet predictions
#'  and additional data used to rerun DoubletFinder with annotation data if needed
#'
#' @export
#'
#' @import Seurat
#' @import DoubletFinder
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' objs <- lapply(objs, runDoubletFinder,
#'   PCs = 10,
#'   sct = FALSE,
#'   cores = 1)
#' }
runDoubletFinder <- function(obj,
                             PCs,
                             sct,
                             cores,
                             ground_truth = NULL,
                             doublet_rate = 0.075,
                             assay = 'RNA'
){
  Seurat::DefaultAssay(obj) <- assay
  sweep.res.list <- DoubletFinder::paramSweep_v3(obj,
                                  PCs = 1:PCs,
                                  sct = sct,
                                  num.cores = cores)
  ## DF internally subsets to 10k cells
  ground_truth <- ground_truth[rownames(sweep.res.list[[1]])]
  if (!is.null(ground_truth)){
    sweep.stats <- DoubletFinder::summarizeSweep(sweep.res.list,
                                  GT = TRUE,
                                  GT.calls = ground_truth)
  } else {
    sweep.stats <- DoubletFinder::summarizeSweep(sweep.res.list, GT = FALSE)
  }

  bcmvn <- DoubletFinder::find.pK(sweep.stats)
  pK <- bcmvn[bcmvn$BCmetric==max(bcmvn$BCmetric),'pK']
  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi <- round(doublet_rate*nrow(obj@meta.data))
  obj <- DoubletFinder::doubletFinder_v3(obj,
                          PCs = 1:PCs,
                          pK = as.numeric(as.character(pK)),
                          nExp = nExp_poi,
                          sct = sct)
  obj@meta.data <- obj@meta.data %>%
    dplyr::rename(DF_pANN = starts_with('pANN')) %>%
    dplyr::rename(DF_classifications = starts_with('DF.classifications'))
  return(obj)
}
