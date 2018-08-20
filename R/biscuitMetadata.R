#' metadata about a Biscuit run using its vcfHeader metadata element
#' 
#' @param x   a bsseq object with a $vcfHeader element
#' 
#' @return    information about the run
#' 
#' @import    VariantAnnotation
#'
#' @export 
biscuitMetadata <- function(bsseq) { 
  if (is.null(metadata(bsseq)$vcfHeader)) {
    stop("metadata(bsseq)$vcfHeader is NULL -- cannot extract anything from it")
  } else {
    meta <- VariantAnnotation::meta(metadata(bsseq)$vcfHeader)
    c("Reference genome"=basename(meta$META["reference","Value"]),
      "Biscuit version"=sub("biscuitV", "", meta$META["source","Value"]),
      "Invocation"=meta$program[,'cmd'])
  }
}
