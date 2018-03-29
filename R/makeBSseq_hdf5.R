#' make an HDF5-backed BSseq object from an (imported) Biscuit BED file
#'
#' @param dt            a data.table (from the BED file) 
#' @param betacols      the beta column names
#' @param covgcols      the coverage column names
#' @param sampleNames   the sample names
#'
#' @return an HDF5-backed BSseq object
#' 
#' @import GenomicRanges
#' @import HDF5Array
#' @import bsseq 
#'
#' @seealso makeBSseq
#'
#' @export 
makeBSseq_hdf5 <- function(dt, betacols, covgcols, sampleNames) {
  hdf5_M <- writeHDF5Array(fixNAs(round(dt[, betacols, with=FALSE]* 
                                        dt[, covgcols, with=FALSE])))
  hdf5_Cov <- writeHDF5Array(fixNAs(dt[, covgcols, with=FALSE]))
  BSseq(gr=makeGRangesFromDataFrame(dt[, c("chr","start","end")]),
        M=hdf5_M, Cov=hdf5_Cov, sampleNames=sampleNames, rmZeroCov=TRUE)
}
