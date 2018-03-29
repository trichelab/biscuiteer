#' a bsseq loader for Biscuit mergecg output (BED-like format, 3 cols/sample)
#' e.g. P01-028-T06.joint.cg.merged.hg19.bed.gz has 3 samples, and thus 12 cols
#'
#' @param filename    the file (compressed or not, doesn't matter) to load
#' @param sampleNames sample names for the bsseq object (if NULL, will create)
#' @param hdf5        make the object HDF5-backed? (FALSE) 
#' 
#' @return            a BSseq object from the bsseq package
#'
#' @import bsseq
#' @import data.table
#' @import HDF5Array
#'
#' @seealso load.biscuit.unmerged
#'
#' @export
load.biscuit.merged <- function(filename, sampleNames=NULL, hdf5=FALSE) {

  params <- checkBiscuitBED(filename, sampleNames, merged=TRUE)
  ncolumns <- 3 + (3 * params$nSamples)
  dropcols <- seq(6, ncolumns, 3)
  merged.dt <- fread(params$input, sep="\t", sep2=",", na.string=".", 
                     drop=dropcols)
  colnames(merged.dt) <- params$colNames
  merged.dt[, "start"] <- merged.dt[, "start"] + 1 # quirk
  message("Loaded data from ", filename, ". Creating bsseq object...")

  if (hdf5) { 
    with(params, makeBSseq_hdf5(merged.dt, betacols, covgcols, sampleNames))
  } else { 
    with(params, makeBSseq(merged.dt, betacols, covgcols, sampleNames))
  } 
}
