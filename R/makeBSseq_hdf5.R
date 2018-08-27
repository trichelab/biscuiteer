#' make an HDF5-backed BSseq object from an (imported) Biscuit BED file
#'
#' @param tbl         a tibble (from read_tsv) or a data.table (from fread())
#' @param params      parameters (from checkBiscuitBED)
#' @param simplify    simplify sample names by dropping .foo.bar.hg19 & similar
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
makeBSseq_hdf5 <- function(tbl, params, simplify=FALSE) { 

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
  hdf5M <- writeHDF5Array(M)
  hdf5Cov <- writeHDF5Array(Cov)
  res <- BSseq(gr=gr, M=hdf5M, Cov=hdf5Cov, pData=params$pData, rmZeroCov=TRUE)
  if (simplify) res <- simplifySampleNames(res)
  return(res)

}
