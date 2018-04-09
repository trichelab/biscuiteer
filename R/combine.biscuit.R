#' a bsseq combiner for Biscuit output (quite possibly sparse-matrix-backed) 
#'
#' @param  x          a BSseq object
#' @param  y          another BSseq object
#' @param  sparse     should the combined object be sparse? (FALSE) 
#' 
#' @return            a combined bsseq::BSseq object
#'
#' @import bsseq
#'
#' @seealso BSseq
#' @seealso load.biscuit.merged
#' @seealso load.biscuit.unmerged
#'
#' @export
combine.biscuit <- function(x, y, sparse=FALSE) {

  message("Preparing to combine biscuit output...")
  if (sparse) {
    stop("Sparse combining is not yet supported")
  } else {
    newgr <- union(granges(x), granges(y))
    newgr <- keepSeqlevels(newgr, 
                           intersect(seqlevels(x), seqlevels(y)),
                           pruning.mode="coarse")
    message("Adding holes...")
    stop("Combination of BSseq objects is not yet finished")
  }

}
