#' Attempt to simplify bsseq sample names (usually derived from VCF columns) 
#' 
#' Tries using the longest common subsequence to figure out what can be dropped.
#' 
#' @param   x     a BSseq (or any SummarizedExperiment-derived) object
#' 
#' @return        the same object but with simplified sample names
#' 
#' @export 
simplifySampleNames <- function(x) { 
  lcs <- function(a, b) {
    paste(qualV::LCS(strsplit(a,"")[[1]], strsplit(b,"")[[1]])$LCS, collapse="")
  }
  sampleNames(x) <- sub(Reduce(lcs, sampleNames(x)), "", sampleNames(x))
  return(x) 
}
