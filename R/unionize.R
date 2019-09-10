#' Combine bsseq objects together without losing information
#'
#' Takes provided bsseq objects, the union of their GRanges, fills out the
#' sites no in the union with 0M/0Cov (or ./. in a sparse Matrix), and returns
#' the even-sparser bsseq holding all of them.
#'
#' All conditions governing cbind(bsseq1, bsseq2, ...) also apply to unionize.
#'
#' @param bs1       A bsseq object (returned unaltered if length(list(...)) = 0)
#' @param ...       One or more bsseq objects to combine with bs1
#' @param parallel  Split the bsseq objects by chrom and parallelize
#'                    (DEFAULT: FALSE)
#' @param onlyChrs  Retain a specific subset of chromosomes - if NULL, keep all
#'                    (DEFAULT: NULL)
#'
#' @return          A larger and more sparse bsseq object
#'
#' @import parallel
#' @import Matrix
#' @import bsseq
#'
#' @seealso unionizeBSseq
#'
#' @examples
#'
#'   tcga_mrg <- system.file("extdata",
#'                           "TCGA_BLCA_A13J_chr11p15_unmerged.bed.gz",
#'                           package = "biscuiteer")
#'   tcga_shf <- system.file("extdata",
#'                           "TCGA_BLCA_A13J_chr11p15_shuffled_unmerged.bed.gz",
#'                           package = "biscuiteer")
#'   tcga_vcf <- system.file("extdata", "TCGA_BLCA_A13J_header_only.vcf.gz",
#'                           package = "biscuiteer")
#'   bisc1 <- read.biscuit(BEDfile = tcga_mrg, VCFfile = tcga_vcf,
#'                         merged = FALSE, genome = "hg38")
#'   bisc2 <- read.biscuit(BEDfile = tcga_shf, VCFfile = tcga_vcf,
#'                         merged = FALSE, genome = "hg38")
#'
#'   comb <- unionize(bisc1, bisc2)
#'
#' @export
#'
unionize <- function(bs1,
                     ...,
                     parallel = FALSE,
                     onlyChrs = NULL) { 

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

  # keep only shared and/or requested chromosomes 
  keptSeqlevels <- intersect(seqlevels(bs1), seqlevels(bs2))
  if (!is.null(onlyChrs)) keptSeqlevels <- intersect(keptSeqlevels, onlyChrs)
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
    
    bigpd <- rbind(colData(bs1), colData(bs2))
    biggrl <- GRangesList(
      bs1=suppressWarnings(GenomicRanges::setdiff(granges(bs1), granges(bs2))),
      bss=suppressWarnings(GenomicRanges::intersect(granges(bs1),granges(bs2))),
      bs2=suppressWarnings(GenomicRanges::setdiff(granges(bs2), granges(bs1)))
    )
    message("  merging methylated read matrices...") 
    M <- as.matrix(unionizeBSseq(bs1, bs2, biggrl, what="M"))
    message("  merging read coverage matrices...") 
    Cov <- as.matrix(unionizeBSseq(bs1, bs2, biggrl, what="Cov"))
    
    # sort is required due to the trick in unionizeBSseq
    sort(BSseq(M, Cov, pData=bigpd, sampleNames=rownames(bigpd), gr=unlist(biggrl)))
  }

}
