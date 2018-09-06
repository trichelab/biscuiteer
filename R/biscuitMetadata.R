#' metadata about a Biscuit run using its vcfHeader metadata element, or its VCF
#' 
#' @param x     a bsseq object with a $vcfHeader element (or NULL)
#' @param VCF   a tabix'ed VCF file (can just be the header) whence x is derived
#' 
#' @return      information about the run
#' 
#' @aliases     getBiscuitMetadata
#' 
#' @import      VariantAnnotation
#'
#' @export 
biscuitMetadata <- function(x=NULL, VCF=NULL) { 

  if (is.null(metadata(x)$vcfHeader) & is.null(VCF)) {
    stop("metadata(x)$vcfHeader is NULL; VCF is NULL: where can it be found?")
  } 
  if (!is.null(metadata(x)$vcfHeader) & !is.null(VCF)) {
    message("You have provided BOTH a BSseq object `x` AND a VCF file `VCF`.") 
    message("If metadata(x)$vcfHeader exists, it will take precedence.")
  } 
  if (is.null(x) | is.null(metadata(x)$vcfHeader) & !is.null(VCF)) {
    vcfHead <- VariantAnnotation::scanVcfHeader(VCF)
    meta <- VariantAnnotation::meta(vcfHead)
  } else {
    meta <- VariantAnnotation::meta(metadata(x)$vcfHeader)
  }
  List("Reference genome"=basename(meta$META["reference","Value"]),
       "Biscuit version"=sub("biscuitV", "", meta$META["source","Value"]),
       "Invocation"=meta$program[,'cmd'])
}

# alias
getBiscuitMetadata <- biscuitMetadata
