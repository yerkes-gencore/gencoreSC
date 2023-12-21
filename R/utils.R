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


#' Subset a Seurat object on multiple variables and values
#'
#' Provide one or more variables to subset on with one or more values to
#'  filter for. This extends the base Seurat subsetting function to allow
#'  multiple subsets to be supplied easily as function arguments. This is
#'  intended for use within other functions, but will still work to return
#'  a subsetted seurat object if you want to use it on it's own, although
#'  Seurat's implementation may be preferred if you aren't constrained to
#'  trying to do multiple subsets in a single function call.
#'
#' @param obj Seurat object
#' @param subset_vars Columns of object metadata to subset on
#' @param subsets Values of `subset_vars` to filter for. If multiple `subset_vars`
#'  are provided, you can specify multiple `subsets` values for each var by
#'  providing a list of character vectors. See examples.
#'
#' @returns A filtered Seurat object
#' @export
#'
#' @examples
#' \dontrun{
#' ## If you only have one value to subset for, the list convention isn't too picky
#' subset_seurat_object(obj = obj,
#'   subset_vars = 'cell_type',
#'   subsets = 'ISG-rich')
#'
#' ## If you want to subset for multiple values of a single subset_var, do this
#' subset_seurat_object(obj = obj,
#'   subset_vars = 'cell_type',
#'   subsets = list(c('ISG-rich', 'Other')))
#'
#' ## If you want to subset for multiple values of multiple subset_vars
#' subset_seurat_object(obj = obj,
#'   subset_vars = c('cell_type', 'timepoint'),
#'   subsets = list(c('ISG-rich', 'Other'), c('pre')))
#' }
subset_seurat_object <- function(obj,
                                 subset_vars,
                                 subsets) {
  if (!missing(subsets) & missing(subset_vars)) {
    stop('Subset_vars not specified for subsetting')
  }
  if (missing(subsets) & !missing(subset_vars)) {
    warning('No subsets specified for subset_vars')
  }
  if (length(subsets) != length(subset_vars)) {
    stop('Error: subsets should be provided as a list of character vectors,
         where the number of vectors equals the number of entries in
         subset_vars. See examples for details')
  }

  if (!missing(subsets) & !missing(subset_vars)) {
    for (i in 1:length(subset_vars)) {
      subset_var <- subset_vars[i]
      if (subset_var %in% colnames(obj@meta.data)) {
        # for (subset in unlist(subsets[[i]])) {
        #   if (subset %in% unique(obj@meta.data[[subset_var]])) {
            obj <- obj[,obj@meta.data[[subset_var]] %in% unlist(subsets[[i]])]
          # } else {
          #   stop('Error: subset value :"', subset,
          #        '" not present as a value of ',
          #        subset_var[i], 'column in object metadata.')
          # }
        # }
      } else {
        stop('Error: subset_var value :"',subset_var,
             '" not present as a column in the object metadata.')
      }
    }
  }
  obj
}
