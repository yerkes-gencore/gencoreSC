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
