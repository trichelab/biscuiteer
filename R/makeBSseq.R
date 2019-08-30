#' Make an in-memory bsseq object from a biscuit BED
#'
#' Beware that any reasonably large BED files may not fit into memory!
#'
#' @param tbl       A tibble (from read_tsv) or a data.table (from fread)
#' @param params    Parameters from checkBiscuitBED
#' @param simplify  Simplify sample names by dropping .foo.bar.hg19? (or
#'                    similar) (DEFAULT: FALSE)
#' @param verbose   Print extra statements? (DEFAULT: FALSE)
#'
#' @return          An in-memory bsseq object
#'
#' @import GenomicRanges
#' @import bsseq
#'
#' @examples
#'
#' @export
#'
makeBSseq <- function(tbl,
                      params,
                      simplify = FALSE,
                      verbose = FALSE) {

  gr <- resize(makeGRangesFromDataFrame(tbl[, 1:3]), 1) 

  # helper fn  
  matMe <- function(x, gr, verbose) {
    if (!is(x, "matrix")) {
      if (verbose) message("Turning a vector into a matrix...")
      x <- as.matrix(x)
    }
    return(x)
  }

  # helper fn  
  fixNames <- function(x, gr, what=c("M","Cov"), verbose=FALSE) {
    if (is.null(rownames(x))) {
      if (verbose) message("Adding rownames...")
      rownames(x) <- as.character(gr)
    }
    colnames(x) <- base::sub("beta", match.arg(what), colnames(x))
    return(x)
  }

  # deal with data.table weirdness 
  if (params$how == "data.table") { 
    betas <- match(params$betaCols, names(tbl))
    covgs <- match(params$covgCols, names(tbl))
    M <- matMe(fixNAs(round(tbl[,..betas]*tbl[,..covgs]),y=0,params$sparse), gr)
    Cov <- matMe(fixNAs(tbl[, ..covgs], y=0, params$sparse), gr)
  } else { 
    M <- with(params, 
              matMe(x=fixNAs(round(tbl[,betaCols]*tbl[,covgCols]), y=0, sparse),
                    gr=gr, verbose=verbose))
    Cov <- with(params, 
                matMe(x=fixNAs(tbl[, covgCols], y=0, sparse), 
                      gr=gr, verbose=verbose))
  }
  Cov <- fixNames(Cov, gr, what="Cov", verbose=verbose)
  M <- fixNames(M, gr, what="M", verbose=verbose)
  colnames(Cov) <- colnames(M) <- params$pData$sampleNames
  if (verbose) message("Creating bsseq object...") 
  res <- BSseq(gr=gr, M=M, Cov=Cov, pData=params$pData,
               rmZeroCov=TRUE, sampleNames=params$pData$sampleNames) 
  if (simplify) res <- simplifySampleNames(res)
  return(res)

}
