#' (e)RRBS settings for dmrseq
#'
#' @param bsseq          A bsseq object
#' @param testCovariate  The pData column to test on
#' @param cutoff         The minimum CpG-wise difference to use (DEFAULT: 0.2)
#' @param bpSpan         Span of smoother AND max gap in DMR CpGs (DEFAULT: 750)
#' @param ...            Other arguments to pass along to dmrseq
#'
#' @return               A GRanges object (same as from dmrseq)
#' 
#' @import dmrseq
#' @import BiocParallel
#'
#' @examples
#'
#'   data(BS.chr21, package="dmrseq")
#'   dat <- BS.chr21
#'
#'   rrbs <- RRBSeq(dat[1:500, ], "Rep", cutoff = 0.05, BPPARAM=BiocParallel::SerialParam())
#'
#' @export
#'
RRBSeq <- function(bsseq,
                   testCovariate,
                   cutoff = 0.2,
                   bpSpan = 750,
                   ...) { 

  dmrseq(filterLoci(bsseq, testCovariate), testCovariate=testCovariate, ..., 
         bpSpan=bpSpan, maxGap=bpSpan, maxGapSmooth=bpSpan*2, cutoff=cutoff)

}

