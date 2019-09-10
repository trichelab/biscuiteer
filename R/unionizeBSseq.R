#' Combine two bsseq objects together
#'
#' Takes two bsseq objects and the overlaps of the union of their rowRanges,
#' fills the sites not in the union with 0M/0Cov (or ./. in a sparse Matrix),
#' then returns a Matrix containing the joint information from both objects.
#'
#' @param bs1     bsseq object number 1
#' @param bs2     bsseq object number 2
#' @param biggrl  A GRangesList with (bs1-only, both, and bs2-only) regions
#' @param what    "M" (methylation) or "Cov" (coverage)
#'
#' @return        A Matrix
#'
#' @import Matrix
#' @import bsseq
#'
#' @seealso unionize
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
#'   grl <- GRangesList(
#'            bs1 = GenomicRanges::setdiff(granges(bisc1), granges(bisc2)),
#'            bss = GenomicRanges::intersect(granges(bisc1), granges(bisc2)),
#'            bs2 = GenomicRanges::setdiff(granges(bisc2), granges(bisc1))
#'          )
#'
#'   M <- unionizeBSseq(bisc1, bisc2, grl, what="M")
#'   #Cov <- unionizeBSseq(bisc1, bisc2, grl, what="Cov")
#'
#' @export
#'
unionizeBSseq <- function(bs1,
                          bs2,
                          biggrl,
                          what = c("M","Cov")) {

  what <- match.arg(what, c("M","Cov"))
  if (!is(bs1, "BSseq")) stop("bs1 must be a BSseq object.")
  if (!is(bs2, "BSseq")) stop("bs2 must be a BSseq object.") 
  if (!identical(seqlevels(bs1), seqlevels(bs2))) stop("Unmatched seqlevels!") 
  if (!all(c("bs1", "bss", "bs2") %in% names(biggrl))) stop("Bogus biggrl!")
  seqlevs <- length(seqlevels(biggrl))

  overall <- sum(sapply(biggrl, length))
  rawdata1 <- Matrix(getCoverage(bs1, regions=c(biggrl$bs1, biggrl$bss),
                                 what="perRegionTotal",
                                 type=what))
  padding1 <- Matrix(0, nrow=length(biggrl$bs2), ncol=ncol(bs1))
  dat1 <- rbind(rawdata1, padding1)
  padding2 <- Matrix(0, nrow=length(biggrl$bs1), ncol=ncol(bs2))
  rawdata2 <- Matrix(getCoverage(bs2, regions=c(biggrl$bss, biggrl$bs2),
                                 what="perRegionTotal", 
                                 type=what))
  dat2 <- rbind(padding2, rawdata2)  
  unionized <- cbind(dat1, dat2)
  return(unionized)

}
