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
#' s <- list(a = "this is a fake suerat_obj", b = seq(1,5))
#' max_min_mean <- function(x) {
#'   result <- c(max = max(s$b), min = min(s$b), mean = mean(s$b))
#'   return(list(s, result))
#' }
#' c(s, result) %<-% max_min_mean(s)
#'
#' # This also works for lapply
#' x <- list(a = 1:10, beta = exp(-3:3), logic = c(TRUE,FALSE,FALSE,TRUE))
#' c(a, beta, logic) %<-% lapply(x, mean)
#' a
#' beta
#' logic
#' }
NULL
