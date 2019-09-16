#' Filter loci with zero coverage
#'
#' Function potentially used to be a part of dmrseq. Included here to avoid
#' dmrseq failing due to any number of reasons related to lack of coverage.
#'
#' The code is adapted from the precheck loop of dmrseq::dmrseq
#'
#' @param bsseq          A bsseq object for filtering
#' @param testCovariate  The name of the pData column dmrseq will test on
#'
#' @return               A bsseq object ready for dmrseq to use
#'
#' @seealso dmrseq
#' @seealso WGBSeq
#' @seealso RRBSeq
#'
#' @importFrom DelayedMatrixStats rowCounts
#' @import bsseq
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
#'   bisc1 <- read.biscuit(BEDfile = shuf_bed, VCFfile = shuf_vcf,
#'                         merged = FALSE)
#'   bisc2 <- read.biscuit(BEDfile = orig_bed, VCFfile = orig_vcf,
#'                         merged = FALSE)
#'
#'   comb <- unionize(bisc1, bisc2)
#'
#'   filt <- filterLoci(comb, "sampleNames")
#'
#' @export
#'
filterLoci <- function(bsseq,
                       testCovariate) {

  filt <- c()
  lev <- unique(pData(bsseq)[[testCovariate]])
  for (l in seq_along(lev)) {
    inLev <- which(pData(bsseq)[[testCovariate]] == lev[l])
    filt <- c(filt, 
              which(DelayedMatrixStats::rowCounts(getCoverage(bsseq[,inLev]),
                                                  value=0) == length(inLev)))
  }
  if (length(filt) > 0) {
    message(length(filt), " loci with 0 coverage in at least 1 condition.")
    retain <- setdiff(seq_len(nrow(bsseq)), filt)
    message("Retaining ", length(retain), " loci.")
    return(bsseq[retain,])
  } else {
    return(bsseq)
  }

} 
