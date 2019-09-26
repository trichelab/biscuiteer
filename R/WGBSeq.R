#' Wrapper for WGBS settings for dmrseq
#'
#' @param bsseq          A bsseq object
#' @param testCovariate  The pData column to test on
#' @param bpSpan         Span of smoother AND 2x max gap in DMR CpGs
#'                         (DEFAULT: 1000)
#' @param ...            Other arguments to pass along to dmrseq
#'
#' @return               A GRanges object (same as from dmrseq)
#'
#' @import dmrseq
#'
#' @examples
#'
#'   data(BS.chr21, package="dmrseq")
#'   dat <- BS.chr21
#'
#'   wgbs <- WGBSeq(dat[1:500, ], "CellType", cutoff = 0.05,
#'                  BPPARAM=SerialParam())
#'
#' @export
#'
WGBSeq <- function(bsseq,
                   testCovariate,
                   bpSpan = 1000,
                   ...) { 
  # FIXME: Why is maxGap = bpSpan / 2 when it is maxGap = bpSpan in RRBSeq.R????

  dmrseq(filterLoci(bsseq, testCovariate), testCovariate=testCovariate, 
         bpSpan=bpSpan, maxGap=(bpSpan/2), maxGapSmooth=bpSpan*2, ...)

}

