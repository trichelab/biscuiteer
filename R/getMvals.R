#' helper function for compartment inference (shrink-by-smoothing)
#' 
#' since we want something resembling M-values for compartment inference,
#' this function grabs suitably deep (>= 3 reads in >=2 sample) measurements
#' and turns them into lightly moderated M-values for compartment calling
#' by Dirichlet smoothing (adding addReads reads to each of M and Covg)
#' 
#' @param x           a BSseq object with methylated and total reads 
#' @param minCov      minimum read coverage for landmarking samples (3)
#' @param minSamp     minimum landmark samples with >= minCov (2)
#' @param k           pseudoreads for smoothing (1)
#' 
#' @return            a matrix of Mvalues, with genomic coordinates as row names
#'
#' @import            bsseq
#'
#' @export
getMvals <- function(x, minCov=3, minSamp=2, k=1) {

  # do any loci in the object have enough read coverage in enough samples? 
  usable <- (DelayedMatrixStats::rowSums2(getCoverage(x) >= minCov) >= minSamp)
  if (!any(usable)) stop("No usable CpG loci ( >= minCov in >= minSamp )!")
    
  # construct a subset of the overall BSseq object with smoothed mvalues 
  getSmoothedMvals(subset(x, usable), k=k, minCov=minCov)

}

# helper fn
getSmoothedMvals <- function(x, k=1, minCov=3, maxFrac=0.5) {

  res <- log2it((getCoverage(x, type="M") + k) / (getCoverage(x) + k + k))
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

# helper fn
log2it <- function(p) log2(p / (1 - p))
