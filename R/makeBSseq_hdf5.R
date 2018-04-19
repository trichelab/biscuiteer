#' make an HDF5-backed BSseq object from an (imported) Biscuit BED file
#'
#' @param tbl           a tibble
#' @param betacols      the beta column names
#' @param covgcols      the coverage column names
#' @param pData         a DataFrame (usually from checkBiscuitBED)
#' @param sparse        make the object Matrix-backed? (TRUE)
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
makeBSseq_hdf5 <- function(tbl, betacols, covgcols, pData, sparse=TRUE) {
  hdf5M <- writeHDF5Array(fixNAs(round(tbl[,betacols]*tbl[,covgcols]),0,sparse))
  hdf5Cov <- writeHDF5Array(fixNAs(tbl[, covgcols], 0, sparse))
  BSseq(gr=makeGRangesFromDataFrame(tbl[, 1:3]), M=hdf5M, Cov=hdf5Cov, 
        pData=pData, rmZeroCov=TRUE)
}
