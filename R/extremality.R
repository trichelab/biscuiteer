#' compute the fraction of a Bernoulli variance achieved by a proportion
#'
#' Note: this works (efficiently) on matrices and DelayedMatrix objects.
#' Also note: since it is posssible for "raw" extremality to be > 1, 
#'            the function does a second pass to correct for this. 
#' 
#' extremal <- c(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0)
#'   milder <- c(0.8, 0.2, 0.1, 0.3, 0.7, 0.1, 0.4, 0.1, 0.1, 0.2)
#' 
#' extremality(extremal, raw=TRUE) 
#' extremality(extremal, raw=FALSE) 
#'
#' FIXME: this breaks when handed a single-row DelayedArray.
#' 
#' @param x    an rectangular object with proportions in it 
#' @param raw  skip the correction pass? (FALSE) 
#' 
#' @return     the extremality of each row (if more than one) of the object 
#' 
#' @importFrom matrixStats rowVars
#' @importFrom matrixStats rowMeans2
#' @importFrom DelayedMatrixStats rowVars
#' @importFrom DelayedMatrixStats rowMeans2
#' 
#' @export 
extremality <- function(x, raw=FALSE) { 

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
