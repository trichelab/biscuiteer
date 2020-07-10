#' Helper function: expanded expit
#'
#' @param x       a vector of values between -Inf and +Inf
#' @param sqz     the amount by which we 'squoze', default is .000001
#'
#' @return        a vector of values between 0 and 1 inclusive
#'
#' @examples
#'
#'   set.seed(1234)
#'   x <- rnorm(n=1000)
#'   summary(x) 
#'
#'   sqz <- 1 / (10**6)
#'   p <- fexpit(x, sqz=sqz)
#'   summary(p)
#'
#'   all( (abs(x - flogit(p)) / x) < sqz )
#'   all( abs(x - flogit(fexpit(x))) < sqz )
#'
#' @export 
fexpit <- function(x, sqz=0.000001) {
  
  midpt <- .5
  squoze <- exp(x)/(1 + exp(x))
  inflate <- 1 / (1 - (sqz * midpt))
  p <- ((squoze - midpt) * inflate) + midpt 
  return(p)
  
}