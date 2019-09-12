#' Compute fraction of a Bernoulli variance
#'
#' Works efficiently on matrices and DelayedMatrix objects. Note that it is
#' possible for "raw" extremality to be greater than 1, so this function does
#' a second pass to correct for this.
#'
#' @param x    A rectangular object with proportions in it
#' @param raw  Skip the correction pass? (DEFAULT: FALSE)
#'
#' @return     The extremality of each row (if more than one) of the object
#'
#' @importFrom matrixStats rowVars
#' @importFrom matrixStats rowMeans2
#' @importFrom DelayedMatrixStats rowVars
#' @importFrom DelayedMatrixStats rowMeans2
#' @importFrom methods is
#'
#' @examples
#'
#'   x <- rnorm(100, mean=0.5, sd=0.15)
#'   x <- matrix(x, nrow=50, ncol=50)
#'
#'   ext <- extremality(x, raw=TRUE)
#'
#' @export
#'
extremality <- function(x,
                        raw = FALSE) { 
# FIXME: this breaks when handed a single-row DelayedArray.
# extremal <- c(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0)
# milder   <- c(0.8, 0.2, 0.1, 0.3, 0.7, 0.1, 0.4, 0.1, 0.1, 0.2)
# 
# extremality(extremal, raw=TRUE) 
# extremality(extremal, raw=FALSE) 

  bernoulliVar <- function(meanx) meanx * (1 - meanx) 
  
  if (is(x, "matrix")) {
    actualVar <- matrixStats::rowVars(x, na.rm=TRUE)
    meanx <- matrixStats::rowMeans2(x, na.rm=TRUE)
  } else if (is(x, "DelayedArray")) {
    actualVar <- DelayedMatrixStats::rowVars(x, na.rm=TRUE)
    meanx <- DelayedMatrixStats::rowMeans2(x, na.rm=TRUE)
  } else { 
    actualVar <- var(x, na.rm=TRUE) 
    meanx <- mean(x, na.rm=TRUE) 
  } 
  
  rawExtr <- actualVar/bernoulliVar(meanx)
  if (is(x, "matrix")) names(rawExtr) <- rownames(x)
  if (is(x, "DelayedArray")) names(rawExtr) <- rownames(x)
  if (raw) return(rawExtr)

  maxExtr <- extremality(round(x), raw=TRUE) 
  adjExtr <- rawExtr / maxExtr

  return(adjExtr) 

}
