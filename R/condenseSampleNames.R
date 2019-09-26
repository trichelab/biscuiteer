#' Simplify sample names for a bsseq object
#'
#' Utility function for extracting sample names from tabix'ed sample columns,
#' assuming a VCF-naming scheme (such as Sample_1.foo, Sample_1.bar or
#' Sample1_foo, Sample1_bar).
#'
#' @param tbx       A TabixFile instance to parse
#' @param stride    How many columns per sample
#' @param trailing  Trailing character to trim (DEFAULT: "\\.$")
#'
#' @return          A character vector of sample names (longest common
#'                    substrings)
#'
#' @importFrom qualV LCS
#'
#' @examples
#'
#'   library(Rsamtools)
#'   orig_bed <- system.file("extdata", "MCF7_Cunha_chr11p15.bed.gz",
#'                           package="biscuiteer")
#'   if (length(headerTabix(orig_bed)$header) > 0) {
#'     condenseSampleNames(orig_bed, 2)
#'   }
#'
#' @export
#'
condenseSampleNames <- function(tbx,
                                stride,
                                trailing = "\\.$") {
  cols <- strsplit(base::sub("^#", "", headerTabix(tbx)$header), "\t")[[1]]
  indexcols <- headerTabix(tbx)$indexColumns
  samplecols <- cols[setdiff(seq_along(cols), indexcols)]
  nSamples <- length(samplecols) / stride
  idxSample <- rep(seq_len(nSamples), each=stride)
  SNs <- sapply(split(samplecols, idxSample), function(x) Reduce(.LCSubstr, x))
  base::gsub(trailing, "", unname(SNs))
}

# Helper function for finding the longest common subsequence
#
# Takes two character vectors - a and b - and returns their longest common
#   subsequence
.LCSubstr <- function(a, b) {
  paste(qualV::LCS(strsplit(a,"")[[1]], strsplit(b,"")[[1]])$LCS, collapse="")
}
