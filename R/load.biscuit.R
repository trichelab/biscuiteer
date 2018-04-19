#' a bsseq loader for Biscuit output (BED-like format, 2 or 3 cols/sample)
#' e.g. P01-028-T06.joint.ch.hg19.bed.gz has 3 samples, and 9 columns total,
#' while P01-028-T06.joint.cg.merged.hg19.bed.gz has 3 samples and 12 columns
#'
#' @param filename    the file (compressed or not, doesn't matter) to load
#' @param sampleNames sample names (if NULL, create; if data.frame, make pData)
#' @param hdf5        make the object HDF5-backed? (FALSE; use in-core storage) 
#' @param sparse      are there a lot of zero-coverage sites? (TRUE, usually)
#' 
#' @return            a bsseq::BSseq object, possibly Matrix- or HDF5-backed
#'
#' @import readr
#' @import bsseq
#'
#' @seealso BSseq
#' @seealso checkBiscuitBED
#'
#' @export
load.biscuit <- function(filename, 
                         sampleNames=NULL, 
                         hdf5=FALSE, 
                         sparse=TRUE) {

  params <- checkBiscuitBED(filename, sampleNames, hdf5=hdf5)
  if (params$passes > 1) { 
    stop("read_tsv_chunked support (for CpH BEDs) is not yet implemented...")
  } else { 
    message("Reading merged CpG input from ", params$tbx$path, "...")
    tbl <- read_tsv(params$tbx$path, na.string=".", comment="#", 
                    col_names=params$colNames, col_types=params$colSpec)
    tbl[, params$colNames[2]] <- tbl[, params$colNames[2]] + 1 # quirk
    message("Loaded ", params$tbx$path, ". Creating bsseq object...")
  }

  if (params$hdf5) { 
    with(params, makeBSseq_hdf5(tbl, betacols, covgcols, pData, sparse))
  } else { 
    with(params, makeBSseq(tbl, betacols, covgcols, pData, sparse))
  } 

}
