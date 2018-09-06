#' Simplify bsseq sample names (usually VCF columns) 
#' 
#' Tries using the longest common subsequence to figure out what can be dropped.
#' 
#' @param   x     a SummarizedExperiment-like object, or a character vector
#' 
#' @return        the same object but with simplified sample names
#' 
#' @export 
simplifySampleNames <- function(x) { 
  lcs <- function(a, b) {
    paste(qualV::LCS(strsplit(a,"")[[1]], strsplit(b,"")[[1]])$LCS, collapse="")
  }
  if (is(x, "SummarizedExperiment")) {
    xx <- sampleNames(x)
  } else { 
    xx <- x 
  } 
  subst <- Reduce(lcs, xx)
  while (!all(grepl(subst, xx)) & nchar(subst) > 2) {
    subst <- substr(subst, 2, nchar(subst))
  }
  if (nchar(subst) < 3) {
    message("Sample names are already simplified. Returning unaltered.")
  } else { 
    if (is(x, "SummarizedExperiment")) {
      sampleNames(x) <- sub(subst, "", sampleNames(x))  
    } else { 
      x <- sub(subst, "", x)
    }
  }
  return(x) 
}
