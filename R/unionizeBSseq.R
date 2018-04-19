#' take two bsseq objects and the overlaps of the union of their rowRanges,
#' pad out un-covered sites with 0M/0Cov (i.e., ./. in a sparse Matrix),
#' and return a Matrix containing the joint information from both. 
#'
#' @param bs1       a BSseq object
#' @param bs2       another BSseq object
#' @param biggrl    a GRangesList with (bs1-only, both, bs2-only) regions
#' @param what      "M" (methylation) or "Cov" (coverage)
#' 
#' @return a DelayedMatrix
#'
#' @import DelayedArray
#' @importFrom Matrix Matrix
#' @importClassesFrom Matrix Matrix 
#' @importMethodsFrom Matrix cbind2 rbind2
#' @import bsseq
#' 
#' @export
unionizeBSseq <- function(bs1, bs2, biggrl, what=c("M","Cov"), parallel=TRUE) {

  what <- match.arg(what, c("M","Cov"))
  if (!is(bs1, "BSseq")) stop("bs1 must be a BSseq object.")
  if (!is(bs2, "BSseq")) stop("bs2 must be a BSseq object.") 
  if (!identical(seqlevels(bs1), seqlevels(bs2))) stop("Unmatched seqlevels!") 
  if (!all(c("bs1", "bss", "bs2") %in% names(biggrl))) stop("Bogus biggrl!")
  seqlevs <- length(seqlevels(biggrl))

  overall <- sum(sapply(biggrl, length))
  rawdata1 <- Matrix(getCoverage(bs1, regions=c(biggrl$bs1, biggrl$bss),
                                 what="perRegionTotal",
                                 type=what), ncol=ncol(bs1))
  padding1 <- Matrix(0, nrow=length(biggrl$bs2), ncol=ncol(bs1))
  dat1 <- rbind(rawdata1, padding1)
  padding2 <- Matrix(0, nrow=length(biggrl$bs1), ncol=ncol(bs2))
  rawdata2 <- Matrix(getCoverage(bs2, regions=c(biggrl$bss, biggrl$bs2),
                                type=what, what="perRegionTotal"), 
                    ncol=ncol(bs2))
  dat2 <- rbind(padding2, rawdata2)  
  unionized <- cbind(dat1, dat2)
  return(unionized)

}
