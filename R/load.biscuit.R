#' a bsseq loader for Biscuit output (BED-like format, 2 or 3 cols/sample)
#' e.g. P01-028-T06.joint.ch.hg19.bed.gz has 3 samples, and 9 columns total,
#' while P01-028-T06.joint.cg.merged.hg19.bed.gz has 3 samples and 12 columns
#'
#' @param filename    the file (compressed or not, doesn't matter) to load
#' @param sampleNames sample names (if NULL, create; if data.frame, make pData)
#' @param hdf5        make the object HDF5-backed? (FALSE; use in-core storage) 
#' @param sparse      are there a lot of zero-coverage sites? (FALSE)
#' 
#' @return            a bsseq::BSseq object, possibly Matrix- or HDF5-backed
#'
#' @import bsseq
#'
#' @seealso BSseq
#' @seealso checkBiscuitBED
#' @seealso load.biscuit.merged
#' @seealso load.biscuit.unmerged
#'
#' @export
load.biscuit <- function(filename, 
                         sampleNames=NULL, 
                         hdf5=FALSE, 
                         sparse=FALSE) {

  params <- checkBiscuitBED(filename, sampleNames)
  if (params$merged) load.biscuit.merged(params) 
  else load.biscuit.unmerged(params)

}
