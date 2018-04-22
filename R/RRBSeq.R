#' wrapper for [e]RRBS settings appropriate to dmrseq
#' 
#' @param bs              a bsseq object
#' @param testCovariate   the pData column to test on 
#' @param bpSpan          span of smoother AND max gap in DMR CpGs (500)
#' @param ...             other arguments to pass along to dmrseq
#'
#' @return                a GRanges object (same as from dmrseq)
#' 
#' @import dmrseq
#'
#' @export
RRBSeq <- function(bs, testCovariate, bpSpan=500, ...) { 

  dmrseq(filterLoci(bs, testCovariate), testCovariate=testCovariate,
         bpSpan=bpSpan, maxGap=bpSpan, maxGapSmooth=bpSpan*2, ...)

}

