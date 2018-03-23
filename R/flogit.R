#' Helper function: squeezed logit
#'
#' @param x       a vector of values in [0, 1] 
#' @param sqz     the amount by which to 'squeeze', default is .000001
#'
#' @return        a vector of values between (-Inf, +Inf) 
#'
#' @import        gtools
#' @export 
flogit <- function(x, sqz=0.000001) {
  x[ which(x < sqz) ] <- sqz 
  x[ which(x > (1 - sqz)) ] <- (1 - sqz)
  logit(x)
}
