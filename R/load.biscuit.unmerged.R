#' a bsseq loader for Biscuit cg/ch output (BED-like format, 2 cols/sample)
#' e.g. P01-028-T06.joint.ch.hg19.bed.gz has 3 samples, and thus 9 columns
#'
#' @param filename    the file (compressed or not, doesn't matter) to load
#' @param sampleNames sample names for the bsseq object (if NULL, will create)
#' @param hdf5        make the object HDF5-backed? (FALSE) 
#' 
#' @return            a BSseq object from the bsseq package
#'
#' @import bsseq
#' @import data.table
#' @import GenomicRanges
#'
#' @concept bsseq
#' @concept BSseq 
#'
#' @seealso load.biscuit
#' @seealso checkBiscuitBED
#'
#' @export
load.biscuit.unmerged <- function(filename, sampleNames=NULL, hdf5=FALSE) {

  params <- checkBiscuitBED(filename, sampleNames, merged=FALSE)
  chh <- ifelse(base::grepl("c(p?)g", filename, ignore=TRUE), "CpG", "CpH")
  message("Reading unmerged ", chh, " input from ", filename, "...")
  unmerged.dt <- fread(params$input, sep="\t", sep2=",", na.string=".") 
  colnames(unmerged.dt) <- params$colNames
  unmerged.dt[, "start"] <- unmerged.dt[, "start"] + 1 # quirk
  message("Loaded data from ", filename, ". Creating bsseq object...")

  if (hdf5) { 
    with(params, makeBSseq_hdf5(unmerged.dt, betacols, covgcols, sampleNames))
  } else { 
    with(params, makeBSseq(unmerged.dt, betacols, covgcols, sampleNames))
  }

}
