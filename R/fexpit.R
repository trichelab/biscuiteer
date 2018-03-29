#' Helper function: expanded expit
#'
#' @param x       a vector of values between -Inf and +Inf
#' @param sqz     the amount by which to 'squeeze', default is .000001
#'
#' @return        a vector of values between 0 and 1 inclusive
#'
#' @import        gtools
#'
#' @export 
fexpit <- function(x, sqz) {
  (((((inv.logit(x) * 2) - 1) / (1 - sqz)) + 1) / 2)
}
