#' getmode function
#'
#' @param v vector
#'
#' @returns returns the mode of the vector
#' @export
#'
#' @examples
#' getmode(c(1,2,2,2,3,4,5,5,6))
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
