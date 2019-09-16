#' Simplify bsseq sample names
#'
#' Tries using the longest common subsequence to figure out what can be
#' dropped. Usually used for VCF columns.
#'
#' @param x  A SummarizedExperiment-derived object, or a character vector
#'
#' @return   The input object, but with simplified sample names
#'
#' @importFrom qualV LCS
#' @importFrom methods is
#'
#' @examples
#'
#'   orig_bed <- system.file("extdata", "MCF7_Cunha_chr11p15.bed.gz",
#'                           package="biscuiteer")
#'   orig_vcf <- system.file("extdata", "MCF7_Cunha_header_only.vcf.gz",
#'                           package="biscuiteer")
#'   bisc <- read.biscuit(BEDfile = orig_bed, VCFfile = orig_vcf,
#'                        merged = FALSE)
#'
#'   bisc <- simplifySampleNames(bisc)
#'
#' @export
#'
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
