#' Pretty much what it says on the tin.  Summarize methylation over regions.
#'
#' @param   bsseq     the BSseq object to summarize 
#' @param   segs      the regions to summarize over 
#' @param   dropNA    whether to drop rows with > half of samples NaN (TRUE)
#' @param   impute    whether to impute NAs (FALSE) 
#'
#' @return  matrix    a matrix of regional methylation fractions
#'
#' @import  DelayedMatrixStats
#' @import  bsseq 
#' @import  impute
#' 
#' @export
summarizeBsSeqOver <- function(bsseq, segs, dropNA=TRUE, impute=FALSE) { 
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
