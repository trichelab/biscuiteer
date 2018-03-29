#' a bsseq loader for Biscuit output (BED-like format, 2 or 3 cols/sample)
#' e.g. P01-028-T06.joint.ch.hg19.bed.gz has 3 samples, and 9 columns total,
#' while P01-028-T06.joint.cg.merged.hg19.bed.gz has 3 samples and 12 columns
#'
#' @param filename    the file (compressed or not, doesn't matter) to load
#' @param sampleNames sample names for the bsseq object (if NULL, will create)
#' @param hdf5        make the object HDF5-backed? (FALSE) 
#' @param merged      is the file a merged CpG file? (if NULL, will guess) 
#' 
#' @return            a bsseq::BSseq object, possibly HDF5-backed
#'
#' @import bsseq
#'
#' @seealso BSseq
#' @seealso checkBiscuitBED
#' @seealso load.biscuit.merged
#' @seealso load.biscuit.unmerged
#'
#' @export
load.biscuit <- function(filename, sampleNames=NULL, hdf5=FALSE, merged=NULL) {

  message("Preparing to load biscuit output from ", filename, "...")
  if (is.null(merged)) merged <- base::grepl("merged", ignore=TRUE, filename)
  if (merged) { 
    load.biscuit.merged(filename, sampleNames, hdf5)
  } else {
    load.biscuit.unmerged(filename, sampleNames, hdf5)
  } 

}
