#' take a bunch of bsseq objects, take the union of their granges(), 
#' pad out un-covered sites with 0M/0Cov (i.e., ./. in a sparse Matrix),
#' and return the now-even-sparser bsseq holding all of them. All 
#' conditions governing cbind(bsseq1, bsseq2) also apply to unionize(). 
#'
#' @param bs1     a bsseq object (return unaltered if length(list(...)) == 0)
#' @param ...     one or more bsseq objects to combine with the first one 
#' @param hdf5    use HDF5 arrays to return an out-of-core object? (FALSE)
#' 
#' @return        a larger and more sparse bsseq object
#' 
#' @import Matrix
#' @import bsseq
#' 
#' @export
unionize <- function(bs1, ..., hdf5=FALSE) { 

  stopifnot(is(bs1, "BSseq"))
  if (length(list(...)) == 0) {
    # return 
    return(bs1)
  } else if (length(list(...)) == 1) {
    # reprocess
    bs2 <- list(...)[[1]]
  } else { 
    # recurse
    bs2 <- Reduce(unionize, list(...))
  }

  # keep only shared
  message("Pruning...") 
  keptSeqLevels <- intersect(seqlevels(bs1), seqlevels(bs2))
  bs1 <- keepSeqlevels(bs1, keptSeqlevels, pruning.mode="coarse")
  bs2 <- keepSeqlevels(bs2, keptSeqlevels, pruning.mode="coarse")

  # this could be made a lot more efficient by not recursing...
  message("Merging colData...") 
  bigpd <- rbind(colData(bs1), colData(bs2))
  message("Merging rowRanges...") 
  biggrl <- GRangesList(
    bs1=suppressWarnings(GenomicRanges::setdiff(granges(bs1), granges(bs2))),
    both=suppressWarnings(GenomicRanges::intersect(granges(bs1), granges(bs2))),
    bs2=suppressWarnings(GenomicRanges::setdiff(granges(bs2), granges(bs1)))
  )
  message("Merging data matrices...") 
  sort(BSseq(M=unionizeBSseq(bs1, bs2, biggrl, what="M"),
             Cov=unionizeBSseq(bs1, bs2, biggrl, what="Cov"),
             pData=bigpd,
             gr=unlist(biggrl)))

}
