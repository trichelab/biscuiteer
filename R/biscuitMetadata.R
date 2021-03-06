#' Biscuit metadata from VCF header
#'
#' Returns metadata from a Biscuit run using either a supplied VCF file or
#' the vcfHeader metadata element from the bsseq object
#'
#' @param bsseq  A bsseq object with a vcfHeader element (DEFAULT: NULL)
#' @param VCF    A tabix'ed VCF file (can just be the header information)
#'                 from which the bsseq vcfHeader element is derived
#'                 (DEFAULT: NULL)
#'
#' @return       Information regarding the Biscuit run
#'
#' @importFrom VariantAnnotation scanVcfHeader meta
#'
#' @aliases getBiscuitMetadata
#'
#' @examples
#'
#'   orig_bed <- system.file("extdata", "MCF7_Cunha_chr11p15.bed.gz",
#'                           package="biscuiteer")
#'   orig_vcf <- system.file("extdata", "MCF7_Cunha_header_only.vcf.gz",
#'                           package="biscuiteer")
#'   bisc <- readBiscuit(BEDfile = orig_bed, VCFfile = orig_vcf,
#'                       merged = FALSE)
#'
#'   meta <- biscuitMetadata(bisc)
#'   
#' @export
#'
biscuitMetadata <- function(bsseq = NULL,
                            VCF = NULL) { 

  if (is.null(metadata(bsseq)$vcfHeader) & is.null(VCF)) {
    stop("metadata(bsseq)$vcfHeader is NULL; ",
         "VCF is NULL: where can it be found?")
  } 
  if (!is.null(metadata(bsseq)$vcfHeader) & !is.null(VCF)) {
    message("You have provided BOTH a BSseq object ",
            "`bsseq` AND a VCF file `VCF`.") 
    message("If metadata(bsseq)$vcfHeader exists, it will take precedence.")
  } 
  if (is.null(bsseq) | is.null(metadata(bsseq)$vcfHeader) & !is.null(VCF)) {
    vcfHead <- VariantAnnotation::scanVcfHeader(VCF)
    meta <- VariantAnnotation::meta(vcfHead)
  } else {
    meta <- VariantAnnotation::meta(metadata(bsseq)$vcfHeader)
  }
  List("Reference genome"=basename(meta$reference[,"Value"]),
       "Biscuit version"=sub("biscuitV", "", meta$source[,"Value"]),
       "Invocation"=meta$program[,'cmd'])
}

#' @describeIn biscuitMetadata Alias for biscuitMetadata
#'
getBiscuitMetadata <- biscuitMetadata
