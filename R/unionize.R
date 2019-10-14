#' Combine bsseq objects together without losing information
#'
#' Wrapper for the combine(bsseq1, ...) method in bsseq
#'
#' Takes provided bsseq objects, the union of their GRanges, fills out the
#' sites not in the union with 0M/0Cov, and returns the even-sparser bsseq
#' holding all of them.
#'
#' @param bs1  A bsseq object
#' @param ...  One or more bsseq objects to combine with bs1
#'
#' @return     A larger and more sparse bsseq object
#'
#' @examples
#'
#'   shuf_bed <- system.file("extdata", "MCF7_Cunha_chr11p15_shuffled.bed.gz",
#'                           package="biscuiteer")
#'   orig_bed <- system.file("extdata", "MCF7_Cunha_chr11p15.bed.gz",
#'                           package="biscuiteer")
#'   shuf_vcf <- system.file("extdata",
#'                           "MCF7_Cunha_shuffled_header_only.vcf.gz",
#'                           package="biscuiteer")
#'   orig_vcf <- system.file("extdata",
#'                           "MCF7_Cunha_header_only.vcf.gz",
#'                           package="biscuiteer")
#'   bisc1 <- readBiscuit(BEDfile = shuf_bed, VCFfile = shuf_vcf,
#'                        merged = FALSE)
#'   bisc2 <- readBiscuit(BEDfile = orig_bed, VCFfile = orig_vcf,
#'                        merged = FALSE)
#'
#'   comb <- unionize(bisc1, bisc2)
#'
#' @export
#'
unionize <- function(bs1,
                     ...) {
  return(combine(bs1,...))
}
