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
  subst <- Reduce(lcs, sampleNames(x))
  while (!all(grepl(subst, sampleNames(x))) & nchar(subst) > 2) {
    subst <- substr(subst, 2, nchar(subst))
  }
  if (nchar(subst) < 3) {
    message("Could not find a long enough common substring to simplify names.")
  } else { 
    sampleNames(x) <- sub(subst, "", sampleNames(x))
  }
  return(x) 
}
