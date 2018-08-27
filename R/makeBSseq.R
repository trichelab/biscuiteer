#' make an in-core BSseq object from a biscuit BED
#'
#' @param tbl         a tibble (from read_tsv) or a data.table (from fread())
#' @param params      parameters (from checkBiscuitBED)
#' @param simplify    simplify sample names by dropping .foo.bar.hg19 & similar
#'
#' @return an in-core BSseq object
#' 
#' @import GenomicRanges
#' @import bsseq 
#'
#' @seealso makeBSseq_HDF5
#'
#' @export 
makeBSseq <- function(tbl, params, simplify=FALSE) {

  gr <- resize(makeGRangesFromDataFrame(tbl[, 1:3]), 1) 
  if (params$how == "data.table") { 
    betas <- match(params$betaCols, names(tbl))
    covgs <- match(params$covgCols, names(tbl))
    M <- fixNAs(round(tbl[, ..betas] * tbl[, ..covgs]), y=0, params$sparse)
    Cov <- fixNAs(tbl[, ..covgs], y=0, params$sparse)
  } else { 
    M <- with(params, fixNAs(round(tbl[,betaCols]*tbl[,covgCols]), y=0, sparse))
    Cov <- with(params, fixNAs(tbl[, covgCols], y=0, sparse)) 
  }
  res <- BSseq(gr=gr, M=M, Cov=Cov, pData=params$pData, rmZeroCov=TRUE) 
  if (simplify) res <- simplifySampleNames(res)
  return(res)

}
