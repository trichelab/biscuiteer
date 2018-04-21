#' filter out loci that have zero coverage for one or more conditions
#' 
#' This function used to be part of dmrseq, or so we recall. In any event, 
#' it is a colossal drag to set up a dmrseq run only to discover that it fails
#' for seemingly mysterious reasons owing to lack of coverage. 
#'
#' The code is adapted from the precheck loop of dmrseq::dmrseq().
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
#' @importFrom DelayedMatrixStats rowCounts
#'
#' @export
filterLoci <- function(bs, testCovariate, ...) { 

  filt <- NULL
  lev <- unique(pData(bs)[[testCovariate]])
  for (l in seq_along(lev)) {
    inLev <- which(pData(bs)[[testCovariate]] == lev[l])
    toDrop <- which(rowCounts(getCoverage(bs[,inLev]),value=0) == length(inLev))
    filt <- c(filt, toDrop)
  }
  if (length(filt) > 0) {
    message(length(filt), " loci with 0 coverage in at least 1 condition.")
    retain <- setdiff(seq_len(nrow(bs)), filt)
    message("Retaining ", length(retain), " loci.")
    return(bs[retain,])
  } else {
    return(bs)
  }

} 
