#' a bsseq loader for Biscuit output (BED-like format, 2 or 3 cols/sample)
#' e.g. P01-028-T06.joint.ch.hg19.bed.gz has 3 samples, and 9 columns total,
#' while P01-028-T06.joint.cg.merged.hg19.bed.gz has 3 samples and 12 columns
#'
#' @param filename    the file (compressed or not, doesn't matter) to load
#' @param sampleNames sample names (if NULL, create; if data.frame, make pData)
#' @param hdf5        make the object HDF5-backed? (FALSE; use in-core storage) 
#' @param merged      is the file a merged CpG file? (if NULL, will guess) 
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
                         merged=NULL, 
                         sparse=FALSE) {

  message("Preparing to load biscuit output from ", filename, "...")
  if (is.null(merged)) merged <- base::grepl("merged", ignore=TRUE, filename)
  if (merged) { 
    load.biscuit.merged(filename=filename,
                        sampleNames=sampleNames, 
                        hdf5=hdf5, 
                        sparse=sparse)
  } else {
    load.biscuit.unmerged(filename=filename,
                          sampleNames=sampleNames,
                          hdf5=hdf5, 
                          sparse=sparse)
  } 

}
