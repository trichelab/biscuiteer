#' make an in-core BSseq object from the results of reading a Biscuit BED file
#'
#' @param dt            a data.table (from the BED file)
#' @param betacols      the beta column names (from checkBiscuitBED) 
#' @param covgcols      the coverage column names (from checkBiscuitBED) 
#' @param pData         a DataFrame (usually from checkBiscuitBED)
#' @param sparse        make the object Matrix-backed? (FALSE)
#'
#' @return an in-core BSseq object
#' 
#' @import GenomicRanges
#' @import data.table
#' @import S4Vectors
#' @import Matrix
#' @import bsseq 
#'
#' @seealso makeBSseq_HDF5
#'
#' @export 
makeBSseq <- function(dt, betacols, covgcols, pData, sparse=FALSE) { 
  BSseq(gr=makeGRangesFromDataFrame(dt[, c("chr","start","end")]),
        M=fixNAs(round(dt[, betacols, with=FALSE]* 
                       dt[, covgcols, with=FALSE]), sparse=sparse),
        Cov=fixNAs(dt[, covgcols, with=FALSE], sparse=sparse), 
        pData=pData, 
        rmZeroCov=TRUE)
}
