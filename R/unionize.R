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
#'   tcga_mrg <- system.file("extdata",
#'                           "TCGA_BLCA_A13J_chr11p15_unmerged.bed.gz",
#'                           package = "biscuiteer")
#'   tcga_shf <- system.file("extdata",
#'                           "TCGA_BLCA_A13J_chr11p15_shuffled_unmerged.bed.gz",
#'                            package = "biscuiteer")
#'   tcga_mvcf <- system.file("extdata",
#'                            "TCGA_BLCA_A13J_header_only.vcf.gz",
#'                            package = "biscuiteer")
#'   tcga_svcf <- system.file("extdata",
#'                            "TCGA_BLCA_A13J_shuffled_header_only.vcf.gz",
#'                            package = "biscuiteer")
#'   bisc1 <- read.biscuit(BEDfile = tcga_mrg, VCFfile = tcga_mvcf,
#'                         merged = FALSE, genome = "hg38")
#'   bisc2 <- read.biscuit(BEDfile = tcga_shf, VCFfile = tcga_svcf,
#'                         merged = FALSE, genome = "hg38")
#'
#'   comb <- unionize(bisc1, bisc2)
#'
#' @export
#'
unionize <- function(bs1,
                     ...) {
  return(combine(bs1,...))
}
