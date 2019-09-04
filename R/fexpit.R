#' Calculate expanded expit function
#'
#' Helper function for calculating expanded expit
#'
#' @param x    A vector of values between -Inf and +Inf
#' @param sqz  The amount by which to 'squeeze' (DEFAULT: 0.000001)
#'
#' @return     A vector of values between 0 and 1 inclusive
#'
#' @importFrom gtools inv.logit
#'
#' @examples
#'
#'   num <- rnorm(100, mean = 0, sd = 100)
#'   exp <- fexpit(num)
#'
#'   num
#'   exp
#'
#' @export
#'
fexpit <- function(x,
                   sqz = 0.000001) {
  tmp <- (gtools::inv.logit(x) * 2) - 1
  tmp <- (tmp / (1 - sqz)) + 1
  tmp <- tmp / 2
  return(tmp)

# Original code
# ((((inv.logit(x) * 2) - 1) / (1 - sqz)) + 1) / 2
}
