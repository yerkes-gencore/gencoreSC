#' Unpack operator
#'
#' See \code{zeallot::\link[zeallot:operator]{\%<-\%}} for details.
#'
#' @name %<-%
#' @rdname unpack
#' @keywords internal
#' @export
#' @importFrom zeallot %<-%
#' @usage x \%<-\% value
#' @param x A name structure, See \code{zeallot::\link[zeallot:operator]{\%<-\%}} for details.
#' @param value A list of values, vector of values, or R object to assign.
#' @return %<-% returns value.
#' This operator is used primarily for its assignment side-effect. %<-% assigns into the environment in which it is evaluated.
#' @examples
#' \dontrun{
#' # Assign list output of a function to two environmental variables
#' s <- list(string = "this is a fake suerat obj", numbers = seq(1,5))
#' max_min_mean <- function(x) {
#'   result <- c(max = max(x$numbers), min = min(x$numbers), mean = mean(x$numbers))
#'   x$new <- result[1]
#'   return(list(x, result))
#' }
#' c(s, result) %<-% max_min_mean(s)
#'
#' # This also works for split seurat objects (sort of)
#' s.split <- list(cap1 = list(string = "this is a fake suerat obj", numbers = seq(1,5)),
#'                 cap2 = list(string = "this is a fake suerat obj", numbers = seq(6,10)))
#'
#' c(c(s.split$cap1, result$cap1),
#'   c(s.split$cap2, result$cap2)) %<-% lapply(s.split, max_min_mean)
#' }
NULL
