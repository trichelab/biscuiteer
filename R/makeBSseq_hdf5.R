#' make an HDF5-backed BSseq object from an (imported) Biscuit BED file
#'
#' @param tbl           a tibble 
#' @param params        parameters (from checkBiscuitBED)
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
makeBSseq_hdf5 <- function(tbl, params) { 
  hdf5M <- 
    writeHDF5Array(fixNAs(round(tbl[,params$betacols]*tbl[,params$covgcols]),
                          y=0, sparse=params$sparse))
  hdf5Cov <- 
    writeHDF5Array(fixNAs(tbl[, params$covgcols], y=0, sparse=params$sparse))
  BSseq(gr=makeGRangesFromDataFrame(tbl[, 1:3]), 
        M=hdf5M, Cov=hdf5Cov, 
        pData=params$pData,
        rmZeroCov=TRUE)
}
