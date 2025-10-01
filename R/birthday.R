#' birthday
#'
#' @param n An integer or vector of integers giving the group size(s).
#'
#' @returns A numeric vector of probabilities, each between 0 and 1,
#' corresponding to the input group sizes.
#' @export
#'
#' @examples
#' #' # Probability for a single class size
#' birthday(23)   # famous case, ~0.507
#'
#' # Probabilities for a range of class sizes
#' birthday(20:25)
#'
#' # For n > 365, probability is always 1
#' birthday(400)
birthday <- function(n) {
  sapply(n, function(k) {
    if (k > 365) return(1)  # pigeonhole principle
    log_p <- lfactorial(365) - lfactorial(365 - k) - k * log(365)
    p_distinct <- exp(log_p)
    1 - p_distinct
  })
}
