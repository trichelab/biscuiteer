#' a bsseq loader for Biscuit mergecg output (BED-like format, 3 cols/sample)
#' e.g. P01-028-T06.joint.cg.merged.hg19.bed.gz has 3 samples, and thus 12 cols
#'
#' @param filename    the file (compressed or not, doesn't matter) to load
#' @param sampleNames sample names (if NULL, create; if data.frame, make pData)
#' @param hdf5        make the object HDF5-backed? (FALSE; use in-core storage) 
#' @param sparse      make the object Matrix-backed? (TRUE; do so if beneficial)
#' 
#' @return            a BSseq object from the bsseq package
#'
#' @import bsseq
#' @import Matrix
#' @import data.table
#' @import HDF5Array
#'
#' @seealso load.biscuit.unmerged
#'
#' @export
load.biscuit.merged <- function(filename, 
                                sampleNames=NULL, 
                                hdf5=FALSE,
                                sparse=TRUE) {

  params <- checkBiscuitBED(filename, sampleNames, merged=TRUE, sparse=sparse)
  ncolumns <- 3 + (3 * params$nSamples)
  dropcols <- seq(6, ncolumns, 3)
  message("Reading merged CpG input from ", filename, "...")
  merged.dt <- fread(params$input, sep="\t", sep2=",", na.string=".", 
                     drop=dropcols, skip=ifelse(hasHeader, 1, 0))
  colnames(merged.dt) <- params$colNames
  merged.dt[, "start"] <- merged.dt[, "start"] + 1 # quirk
  message("Loaded data from ", filename, ". Creating bsseq object...")

  if (hdf5) { 
    with(params, 
         makeBSseq_hdf5(merged.dt, 
                        betacols, 
                        covgcols, 
                        pData, 
                        sparse))
  } else { 
    with(params, 
         makeBSseq(merged.dt, 
                   betacols, 
                   covgcols, 
                   pData,
                   sparse))
  } 
}
