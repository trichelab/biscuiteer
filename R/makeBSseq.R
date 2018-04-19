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
  BSseq(gr=makeGRangesFromDataFrame(tbl[, 1:3]),
        M=fixNAs(round(tbl[,params$betacols]*tbl[,params$covgcols]), 
                 y=0, params$sparse),
        Cov=fixNAs(tbl[, params$covgcols], y=0, params$sparse), 
        pData=params$pData, 
        rmZeroCov=TRUE)
}
