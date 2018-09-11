#' BSseq class methods (VCF-centric) added by biscuiteer (centralized manpage).
#'
#'
#' biscuiteer adds VariantAnnotation methods to BSseq objects with VCF headers:
#' `samples`,`header`,`meta`,`fixed`,`info`,`geno`
#' 
#' Due to inherited method signatures, the argument (singular) to the method may
#' be named `x` or it may be named `object`. Either way, it is a BSseq object.
#'
#'
#' These add to the existing methods defined in package bsseq for class BSseq:
#' `[`,`length`,`sampleNames`,`sampleNames<-`,`pData`,`pData<-`,`show`,`combine`
#' 
#' Those add to the methods BSseq inherits from SummarizedExperiment, such as:
#' `colData`,`rowRanges`,`metadata`,`subset`,`subsetByOverlaps`,`isDisjoint`,&c.
#'
#'
#' Most of the biscuiteer methods operate on the VCF header, which read.biscuit
#' likes to stuff into the `metadata` slot of BSseq objects it produces. Some 
#' may be handy for populating a BSseq object with QC stats, or querying those. 
#' 
#' @param object a BSseq object, preferably with !is.null(metadata(x)$vcfHeader)
#' @param x      a BSseq object, preferably with !is.null(metadata(x)$vcfHeader)
#'
#' @return    depends on the method; usually a List-like object of some sort 
#' 
#' @name      biscuiteer-methods
#' @aliases   BSseq-methods
#' @aliases   coverage
#' @aliases   header
#' @aliases   reference
#' 
#' @seealso   RangedSummarizedExperiment
#' @seealso   VCFHeader-class
#' @seealso   BSseq-class
#' @seealso   BSseq
NULL


#' @rdname biscuiteer-methods
#' @export
setMethod("samples", "BSseq", function(object) samples(.mdvh(x)))


#' @rdname biscuiteer-methods
#' @export
setMethod("header", "BSseq", function(x) header(.mdvh(x)))


#' @rdname biscuiteer-methods
#' @export
setMethod("meta", "BSseq", function(x) meta(.mdvh(x)))


#' @rdname biscuiteer-methods
#' @export
setMethod("fixed", "BSseq", function(x) fixed(.mdvh(x)))


#' @rdname biscuiteer-methods
#' @export
setMethod("info", "BSseq", function(x) info(.mdvh(x)))


#' @rdname biscuiteer-methods
#' @export
setMethod("geno", "BSseq", function(x) geno(.mdvh(x)))


# helper fn
.mdvh <- function(x) { 
  mdvh <- metadata(x)$vcfHeader
  if (is.null(mdvh)) {
    message("metadata(x)$vcfHeader is empty (where `x` is your BSseq object).")
    return(VariantAnnotation::VCFHeader())
  } else {
    return(mdvh)
  }
} 
