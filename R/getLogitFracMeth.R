#' helper function for compartment inference (shrink-by-smoothing logit frac5mC)
#' 
#' We want something with nominally Gaussian error for compartment inference, so
#' this function grabs suitable (default >= 3 reads in >=2 sample) measurements
#' and turns them into lightly moderated, logit-transformed methylated-fraction
#' estimates (also known, unfortunately, as M-values) for compartment calling,
#' by performing Dirichlet smoothing (adding `k` reads to M and U support).
#' 
#' @param x           a BSseq object with methylated and total reads 
#' @param minCov      minimum read coverage for landmarking samples (3)
#' @param minSamp     minimum landmark samples with >= minCov (2)
#' @param k           pseudoreads for smoothing (0.1)
#' 
#' @return            smoothed logit(M/Cov) matrix with coordinates as row names
#'
#' @aliases           getMvals
#' 
#' @import            gtools
#' @import            bsseq
#'
#' @export
getLogitFracMeth <- function(x, minCov=3, minSamp=2, k=0.1) {

  # do any loci in the object have enough read coverage in enough samples? 
  usable <- (DelayedMatrixStats::rowSums2(getCoverage(x) >= minCov) >= minSamp)
  if (!any(usable)) stop("No usable CpG loci ( >= minCov in >= minSamp )!")
    
  # construct a subset of the overall BSseq object with smoothed mvalues 
  getSmoothedLogitFrac(subset(x, usable), k=k, minCov=minCov)

}

# helper fn
getSmoothedLogitFrac <- function(x, k=0.1, minCov=3, maxFrac=0.5) {

  res <- logit((getCoverage(x, type="M") + k) / (getCoverage(x) + k + k))
  rownames(res) <- as.character(granges(x))
  makeNA <- getCoverage(x) < minCov 
  maxPct <- paste0(100 * maxFrac, "%")
  tooManyNAs <- (DelayedMatrixStats::colSums2(makeNA)/nrow(x)) > maxFrac
  if (any(tooManyNAs)) {
    message(paste(colnames(x)[tooManyNAs],collapse=", ")," are >",maxPct," NA!")
  }
  res[ makeNA ] <- NA
  return(res)

}

# alias 
getMvals <- getLogitFracMeth
