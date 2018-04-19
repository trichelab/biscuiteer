#' take a bunch of bsseq objects, take the union of their granges(), 
#' pad out un-covered sites with 0M/0Cov (i.e., ./. in a sparse Matrix),
#' and return the now-even-sparser bsseq holding all of them. All 
#' conditions governing cbind(bsseq1, bsseq2) also apply to unionize(). 
#'
#' @param bs1       a bsseq object (return unaltered if length(list(...)) == 0)
#' @param ...       one or more bsseq objects to combine with the first one 
#' @param parallel  split the bsseq objects by chrom and parallelize? (FALSE) 
#' 
#' @return        a larger and more sparse bsseq object
#'
#' @import parallel 
#' @importFrom Matrix Matrix
#' @importClassesFrom Matrix Matrix 
#' @importMethodsFrom Matrix cbind2 rbind2
#' @import methods 
#' @import bsseq
#' 
#' @export
unionize <- function(bs1, ..., parallel=FALSE) { 

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
  keptSeqlevels <- intersect(seqlevels(bs1), seqlevels(bs2))
  bs1 <- keepSeqlevels(bs1, keptSeqlevels, pruning.mode="coarse")
  bs2 <- keepSeqlevels(bs2, keptSeqlevels, pruning.mode="coarse")

  # worker function 
  unionizeChrom <- function(sub1, sub2) {
    chrom <- unique(unique(seqnames(sub1)), unique(seqnames(sub2)))
    if (length(chrom) > 1) stop("unionizeChrom should never see 2+ chroms")
    message("Unionizing ", chrom, "...")
    unionize(sub1, sub2, parallel=FALSE)
  }

  chroms <- unique(unique(seqnames(bs1)), unique(seqnames(bs2)))
  if (length(chroms) > 1) {
    if (parallel) {
      sort(do.call(rbind, 
                   mcmapply(unionizeChrom, 
                            split(bs1, seqnames(bs1)),
                            split(bs2, seqnames(bs2)))))
    } else { 
      sort(do.call(rbind, 
                   mapply(unionizeChrom, 
                          split(bs1, seqnames(bs1)),
                          split(bs2, seqnames(bs2)))))
    }
  } else { 
    
    message("Merging colData...") 
    bigpd <- rbind(colData(bs1), colData(bs2))

    message("Merging rowRanges...") 
    biggrl <- GRangesList(
      bs1=suppressWarnings(GenomicRanges::setdiff(granges(bs1), granges(bs2))),
      bss=suppressWarnings(GenomicRanges::intersect(granges(bs1),granges(bs2))),
      bs2=suppressWarnings(GenomicRanges::setdiff(granges(bs2), granges(bs1)))
    )

    message("Merging methylated read matrices...") 
    M <- unionizeBSseq(bs1, bs2, biggrl, what="M", parallel=parallel)
    
    message("Merging read coverage matrices...") 
    Cov <- unionizeBSseq(bs1, bs2, biggrl, what="Cov", parallel=parallel)
    
    # sort is required due to the trick in unionizeBSseq
    sort(BSseq(M, Cov, pData=bigpd, gr=unlist(biggrl)))
  }

}
