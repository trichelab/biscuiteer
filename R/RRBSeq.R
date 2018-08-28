#' [e]RRBS settings appropriate to dmrseq
#' 
#' @param bs              a bsseq object
#' @param testCovariate   the pData column to test on 
#' @param cutoff          the minimum CpG-wise difference to use (0.2)
#' @param bpSpan          span of smoother AND max gap in DMR CpGs (750)
#' @param ...             other arguments to pass along to dmrseq
#'
#' @return                a GRanges object (same as from dmrseq)
#' 
#' @import dmrseq
#'
#' @export
RRBSeq <- function(bs, testCovariate, cutoff=0.2, bpSpan=750, ...) { 

  dmrseq(filterLoci(bs, testCovariate), testCovariate=testCovariate, ..., 
         bpSpan=bpSpan, maxGap=bpSpan, maxGapSmooth=bpSpan*2, cutoff=cutoff)

}

