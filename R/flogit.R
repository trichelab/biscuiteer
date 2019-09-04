#' Calculate squeezed logit function
#'
#' Helper function for calculating squeezed logit
#'
#' @param x    A vector of values between 0 and 1 inclusive
#' @param sqz  The amount by which to 'squeeze' (DEFAULT: 0.000001)
#'
#' @return     A vector of values between -Inf and +Inf
#'
#' @importFrom gtools logit
#'
#' @examples
#'
#'   num <- runif(100, min = 0, max = 1)
#'   log <- flogit(num)
#'
#'   num
#'   log
#'
#' @export 
#'
flogit <- function(x,
                   sqz = 0.000001) {
  x[ which(x < sqz) ] <- sqz 
  x[ which(x > (1 - sqz)) ] <- (1 - sqz)
  gtools::logit(x)
}
