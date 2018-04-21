#' filter out loci that have zero coverage for one or more conditions
#' 
#' This function used to be part of dmrseq, or so we recall. In any event, 
#' it is a colossal drag to set up a dmrseq run only to discover that it fails
#' for seemingly mysterious reasons owing to lack of coverage. 
#'
#' The code is adapted directly from the precheck loop of dmrseq::dmrseq().
#'
#' @param bs              a bsseq object for filtration
#' @param testCovariate   the name of the pData column dmrseq will test on 
#' @param ...             not currently used; may pass to dmrseq in the future
#'
#' @return                a bsseq object ready for dmrseq to use
#'
#' @seealso dmrseq
#' @seealso WGBSeq
#' @seealso RRBSeq
#'
#' @import bsseq
#' @importFrom DelayedMatrixStats rowSums2
#'
#' @export
filterLoci <- function(bs, testCovariate, ...) { 

  filter <- NULL
  lev <- unique(pData(bs)[[testCovariate]])
  for (l in seq_along(lev)) {
    inLev <- pData(bs)[[testCovariate]] == lev[l]
    toDrop <- 1 * DelayedMatrixStats::rowSums2(getCoverage(bs[,inLev]) == 0)
    filter <- rbind(filter, toDrop)
  }
  filter <- which(apply(filter, 2, max) > 0)
  if (length(filter) > 0) {
    message(length(filter), " loci with 0 coverage in at least 1 condition.")
    retain <- setdiff(seq_len(nrow(bs)), filter)
    message("Retaining ", length(retain), " loci.")
    return(bs[retain,])
  } else {
    return(bs)
  }

} 
