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

  gr <- resize(makeGRangesFromDataFrame(tbl[, 1:3]), 1)
  if (params$how == "data.table") { 
    betas <- match(params$betacols, names(tbl))
    covgs <- match(params$covgcols, names(tbl))
    M <- fixNAs(round(tbl[, ..betas] * tbl[, ..covgs]), y=0, params$sparse)
    Cov <- fixNAs(tbl[, ..covgs], y=0, params$sparse)
  } else { 
    M <- with(params, fixNAs(round(tbl[,betacols]*tbl[,covgcols]), y=0, sparse))
    Cov <- with(params, fixNAs(tbl[, covgcols], y=0, sparse)) 
  } 
  hdf5M <- writeHDF5Array(M)
  hdf5Cov <- writeHDF5Array(Cov)
  BSseq(gr=gr, M=hdf5M, Cov=hdf5Cob, pData=params$pData, rmZeroCov=TRUE)

}
