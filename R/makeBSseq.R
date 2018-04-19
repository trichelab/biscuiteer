#' make an in-core BSseq object from the results of reading a Biscuit BED file
#'
#' @param tbl           a tibble (from read_tsv)
#' @param betacols      the beta column names (from checkBiscuitBED) 
#' @param covgcols      the coverage column names (from checkBiscuitBED) 
#' @param pData         a DataFrame (usually from checkBiscuitBED)
#' @param sparse        make the object Matrix-backed? (TRUE)
#'
#' @return an in-core BSseq object
#' 
#' @import GenomicRanges
#' @import bsseq 
#'
#' @seealso makeBSseq_HDF5
#'
#' @export 
makeBSseq <- function(tbl, betacols, covgcols, pData, sparse=TRUE) { 
  BSseq(gr=makeGRangesFromDataFrame(tbl[, 1:3]),
        M=fixNAs(round(tbl[, betacols] * tbl[, covgcols]), y=0, sparse=sparse),
        Cov=fixNAs(tbl[, covgcols], y=0, sparse=sparse), 
        pData=pData, rmZeroCov=TRUE)
}
