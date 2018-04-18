#' utility function for extracting samplenames from Tabix'ed sample columns
#' assumes VCF-style naming, such as with Sample_1.foo, Sample_1.bar, 
#' or Sample1_foo,Sample1_bar, or really anything along these lines
#'
#' @param     tbx       a TabixFile instance to parse
#' @param     stride    how many columns per sample 
#' @param     trailing  trailing character to trim ("\\.$")
#'
#' @return    a character vector of sample names (longest common substrings)
#'
#' @import    qualV
#'
#' @export
condenseSampleNames <- function(tbx, stride, trailing="\\.$") {
  cols <- strsplit(base::sub("^#", "", headerTabix(tbx)$header), "\t")[[1]]
  indexcols <- headerTabix(tbx)$indexColumns
  samplecols <- cols[setdiff(seq_along(cols), indexcols)]
  nSamples <- length(samplecols) / stride
  idxSample <- rep(seq_len(nSamples), each=stride)
  SNs <- sapply(split(samplecols, idxSample), function(x) Reduce(.LCSubstr, x))
  base::gsub(trailing, "", unname(SNs))
}

.LCSubstr <- function(a, b) {
  paste(qualV::LCS(strsplit(a,"")[[1]], strsplit(b,"")[[1]])$LCS, collapse="")
}
