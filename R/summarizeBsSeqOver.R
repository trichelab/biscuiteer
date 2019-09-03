#' Summarize methylation over provided regions
#'
#' Used for bsseq objects. Mostly a local wrapp for getMeth.
#'
#' @param bsseq   The bsseq object to summarize
#' @param segs    Regions to summarize over (GRanges object, no GRangesList yet)
#' @param dropNA  Whether to drop rows if more than half of samples are NA
#'                  (DEFAULT: FALSE)
#' @param impute  Whether to impute NAs/NaNs (DEFAULT: FALSE)
#'
#' @return        A matrix of regional methylation fractions
#'
#' @importFrom DelayedMatrixStats rowSums2
#' @import impute
#' @import bsseq
#'
#' @examples
#'
#' @export
#'
summarizeBsSeqOver <- function(bsseq,
                               segs,
                               dropNA = FALSE,
                               impute = FALSE) { 
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
