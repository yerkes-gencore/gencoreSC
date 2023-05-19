#' Extract SoupX contamination correction values for a gene
#'
#' Pull unmodified transcript counts, soupX corrected counts, portion of
#'  cells set to 0 expression, and soup profile counts for a single gene
#'  from a modified SoupChannel produced by `runSoupX()`. This can be used
#'  to see if a gene of interest from a DGE or similar analysis had strong
#'  presence in the background contamination and may be influenced by SoupX
#'  correction.
#'
#' @param gene A gene symbol
#' @param sc SoupChannel object produced by `runSoupX()`
#'
#' @return A list of named numeric values
#' @export
#'
#' @examples
#' \dontrun{
#' sc <- runSoupX(tod = tod,
#'   toc = toc,
#'   clus = clusters)
#' extractGeneCorrectionsSoupX('ZFPM2', sc)
#' }
extractGeneCorrectionsSoupX <- function(gene, sc){
  cntSoggy = sum(sc$toc[gene,] > 0)
  cntStrained = sum(sc$adjusted_counts[gene,] > 0)
  ratio = (cntSoggy - cntStrained) / cntSoggy

  list(
    'Base capture transcript count' = sum(sc$toc[gene,]),
    'SoupX corrected capture count' = sum(sc$adjusted_counts[gene,]),
    'Portion of cells set to 0'     = round(ratio,3),
    'Soup profile transcript count' = sc$soupProfile[gene,'counts']
  )
}
