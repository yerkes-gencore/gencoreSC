#' Write differential expression test results to excel sheet
#'
#' Take a list of differential expression results tables returned from `Seurat::FindMarkers()`
#'  or a similar function and write them to an excel worksheet, with one sheet for each result.
#'
#' @param results List of outputs from `Seurat::FindMarkers()` or a similar function
#' @param sheet_names List of names for sheets, defaults to names of `results`
#' @param output_name Name of file
#' @param outdir Location of file
#'
#' @return `NULL`
#' @export
#'
#' @import openxlsx
#'
#' @examples
#' \dontrun{
#' dge_results$NK_adaptive.post.1D3_vs_cnt <- FindMarkers(object = cd8_obj,
#'                                                        group.by = 'condition.timepoint',
#'                                                        ident.1 = 'post-ATI.1D3', ident.2 = 'post-ATI.control',
#'                                                        subset.ident = 'NK (adaptive)',
#'                                                        features = dge_genes,
#'                                                        logfc.threshold = 0, min.pct = 0,
#'                                                        test.use = 'MAST', assay = 'RNA')
#' dge_results$NK_adaptive.pre.1D3_vs_cnt <- FindMarkers(object = cd8_obj,
#'                                                       group.by = 'condition.timepoint',
#'                                                       ident.1 = 'pre-ATI.1D3', ident.2 = 'pre-ATI.control',
#'                                                       subset.ident = 'NK (adaptive)',
#'                                                       features = dge_genes,
#'                                                       logfc.threshold = 0, min.pct = 0,
#'                                                       test.use = 'MAST', assay = 'RNA')
#' writeDGEResults(dge_results)
#' }
writeDGEResults <- function(results,
                            sheet_names = names(results),
                            output_name = paste0("differentially_expressed_genes.xlsx"),
                            outdir = here::here('outputs')){
  outfile <- file.path(outdir, output_name)
  message(paste0('Writing results to ', outfile))
  wb <- openxlsx::createWorkbook('ENPRC Gencore')
  mapply(FUN = .addWorksheet_DGEres,
         result = results,
         sheet_name = sheet_names,
         MoreArgs = list(wb = wb)
  )
  openxlsx::saveWorkbook(wb, outfile, overwrite = TRUE)
}

.addWorksheet_DGEres <- function(wb,
                                 result,
                                 sheet_name){
  if (nchar(sheet_name) > 31) {
    sheet_name <- substr(sheet_name, 1, 31)
  }
  result <- result[order(result$p_val_adj),]
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb,
                      sheet = sheet_name,
                      x = as.data.frame(result),
                      rowNames = TRUE)
}
