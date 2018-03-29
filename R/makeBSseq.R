#' make an in-core BSseq object from the results of reading a Biscuit BED file
#'
#' @param dt            a data.table (from the BED file)
#' @param betacols      the beta column names
#' @param covgcols      the coverage column names
#' @param sampleNames   the sample names
#'
#' @return an in-core BSseq object
#' 
#' @import GenomicRanges
#' @import bsseq 
#'
#' @seealso makeBSseq_HDF5
#'
#' @export 
makeBSseq <- function(dt, betacols, covgcols, sampleNames) { 
  BSseq(gr=makeGRangesFromDataFrame(dt[, c("chr","start","end")]),
        M=fixNAs(round(dt[, betacols, with=FALSE]* 
                       dt[, covgcols, with=FALSE])),
        Cov=fixNAs(dt[, covgcols, with=FALSE]), 
        sampleNames=sampleNames, 
        rmZeroCov=TRUE)
}
