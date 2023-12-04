#' Get Seurat object `Ident` name
#'
#' Gets the metadata column name of the current Ident of the Seurat object
#'
#' Will error if the `Ident` isn't derived from a metadata column. Not sure why that would ever happen though...
#'
#' @param obj Seurat object
#'
#' @return Name of metadata column used as the active.ident.
#'
#' @examples
#' \dontrun{
#' # Get name of Ident
#' IdentName(obj)
#' }
#'
#' @export
IdentName <- function(obj) {
  Ident.name <- lapply(obj@meta.data, FUN = function(md, aI) {
    identical(md, unname(aI))
  }, aI = obj@active.ident
  ) %>% unname() %>% unlist() %>% which() %>% names(obj@meta.data)[.]

  if (identical(Ident.name, character(0))) {
    errorCondition("Ident is not in the metadata.")
  }

  return(Ident.name)
}

#' Read an excel workbook with sheets to a list of tables
#'
#' Read an excel workbook with sheets to a list of tables
#'
#' @param filename File in format .xls or .xlsx
#' @param tibble BOOL, Return the data as a tibble instead of a data.frame
#' @param \dots Additional arguments passed to `openxlsx::read.xlsx()`
#'
#' @returns A list of data.frames (or tibbles)
#'
#' @import openxlsx
#' @export
#'
#' @examples
#' \dontrun{
#'   dge_results <- read_excel_allsheets(here('outputs/dge.xlsx'))
#' }
read_excel_allsheets <- function(filename,
                                 tibble = FALSE,
                                 ...) {
  # https://stackoverflow.com/a/12945838/15664425
  sheets <- openxlsx::getSheetNames(filename)
  x <- lapply(sheets, function(X) openxlsx::read.xlsx(filename, sheet = X, ...))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}
