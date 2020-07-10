#' Helper function: squeezed logit
#'
#' @param p       a vector of values between 0 and 1 inclusive
#' @param sqz     the amount by which to 'squeeze', default is .000001
#'
#' @return        a vector of values between -Inf and +Inf
#'
#' @examples
#'
#'   set.seed(1234)
#'   p <- runif(n=1000)
#'   summary(p) 
#' 
#'   sqz <- 1 / (10**6)
#'   x <- flogit(p, sqz=sqz)
#'   summary(x) 
#'
#'   all( abs(p - fexpit(x, sqz=sqz)) < sqz )
#'   all( abs(p - fexpit(flogit(p, sqz=sqz), sqz=sqz)) < sqz ) 
#'
#' @export 
flogit <- function(p, sqz=0.000001) { 
  
  midpt <- 0.5
  deflate <- 1 - (sqz * midpt)
  if (any(p > 1 | p < 0, na.rm = TRUE)) stop("Values of p outside (0,1) detected.")
  squoze <- ((p - midpt) * deflate) + midpt
  return( log( squoze / (1 - squoze)) )
  
}