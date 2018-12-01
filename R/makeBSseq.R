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

  # helper fn  
  matMe <- function(x, gr) {
    if (!is(x, "matrix")) return(as.matrix(x)) 
    else return(x)
  }

  # deal with data.table weirdness 
  if (params$how == "data.table") { 
    betas <- match(params$betaCols, names(tbl))
    covgs <- match(params$covgCols, names(tbl))
    M <- matMe(fixNAs(round(tbl[,..betas]*tbl[,..covgs]),y=0,params$sparse), gr)
    Cov <- matMe(fixNAs(tbl[, ..covgs], y=0, params$sparse), gr)
  } else { 
    M <- with(params, 
              matMe(fixNAs(round(tbl[,betaCols]*tbl[,covgCols]),y=0,sparse),gr))
    Cov <- with(params, 
                matMe(fixNAs(tbl[, covgCols], y=0, sparse), gr))
  }
  colnames(M) <- base::sub("beta", "M", colnames(M))
  if (is.null(rownames(M))) rownames(M) <- as.character(gr)
  colnames(Cov) <- base::sub("beta", "Cov", colnames(Cov))
  if (is.null(rownames(Cov))) rownames(Cov) <- as.character(gr)
  res <- BSseq(gr=gr, M=M, Cov=Cov, pData=params$pData, rmZeroCov=TRUE) 
  if (simplify) res <- simplifySampleNames(res)
  return(res)

}
