#' make an in-core BSseq object from the results of reading a Biscuit BED file
#'
#' @param tbl           a tibble (from read_tsv)
#' @param params        parameters (from checkBiscuitBED)
#'
#' @return an in-core BSseq object
#' 
#' @import GenomicRanges
#' @import bsseq 
#'
#' @seealso makeBSseq_HDF5
#'
#' @export 
makeBSseq <- function(tbl, params) {

  gr <- resize(makeGRangesFromDataFrame(tbl[, 1:3]), 1) 
  if (params$how == "data.table") { 
    M <- with(params, 
              fixNAs(round(tbl[, ..betacols] * tbl[, ..covgcols]), y=0, sparse))
    Cov <- with(params,
                fixNAs(tbl[, ..covgcols], y=0, sparse)) 
  } else { 
    M <- with(params, fixNAs(round(tbl[,betacols]*tbl[,covgcols]), y=0, sparse))
    Cov <- with(params, fixNAs(tbl[, covgcols], y=0, sparse)) 
  } 
  BSseq(gr=gr, M=M, Cov=Cov, pData=params$pData, rmZeroCov=TRUE) 

}
