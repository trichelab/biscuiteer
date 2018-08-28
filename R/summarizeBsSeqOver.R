#' Pretty much what it says on the tin.
#' 
#' Summarize methylation over regions. Mostly a wrapper for getMeth. 
#'
#' @param   bsseq     the BSseq object to summarize 
#' @param   segs      regions to summarize over (a GRanges; no GRangesLists yet)
#' @param   dropNA    whether to drop rows with > half of samples NaN (FALSE)
#' @param   impute    whether to impute NAs/NaNs (FALSE) 
#'
#' @return  matrix    a matrix of regional methylation fractions
#'
#' @importFrom  DelayedMatrixStats rowSums2
#' @import  impute
#' @import  bsseq 
#' 
#' @export
summarizeBsSeqOver <- function(bsseq, segs, dropNA=FALSE, impute=FALSE) { 
  segs <- subsetByOverlaps(segs, bsseq)
  res <- bsseq::getMeth(bsseq, regions=segs, what="perRegion", type="raw")
  rownames(res) <- as.character(segs)
  if (dropNA) {
    res <- res[DelayedMatrixStats::rowSums2(is.na(res)) < (ncol(res)/2),]
  }
  if (impute && any(is.nan(res))) {
    res <- DelayedArray(fexpit(impute.knn(flogit(as.matrix(res)))$data))
  }
  return(res)
} 
