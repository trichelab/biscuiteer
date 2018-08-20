#' metadata about a Biscuit run using its vcfHeader metadata element
#' 
#' @param x   a bsseq object with a $vcfHeader element
#' 
#' @return    information about the run
#' 
#' @import    VariantAnnotation
#'
#' @export 
biscuitMetadata <- function(x) { 
  if (is.null(metadata(x)$vcfHeader)) {
    stop("metadata(x)$vcfHeader is NULL -- cannot extract anything from it")
  } else {
    meta <- VariantAnnotation::meta(metadata(x)$vcfHeader)
    res <- c("Reference genome"=basename(meta$META["reference","Value"]),
             "Biscuit version"=sub("biscuitV", "", meta$META["source","Value"]),
             "Invocation"=meta$program[,'cmd'])
    dat <- DataFrame(Value=res)
    rownames(dat) <- names(res)
    return(dat)
  }
}
