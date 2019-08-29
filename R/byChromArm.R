#' A simple parallization step
#'
#' This function splits an object by chromosome arm, which tends to make
#' parallelization much easier, as cross-arm dependencies are unusual.
#' Therefore, the larger chromosomes can be split across processes or machines
#' without worrying much about data starvation for processes on smaller
#' chromosomes.
#'
#' @param x     Any object with a GRanges in it: bsseq, SummarizedExperiment...
#' @param arms  Another GRanges, but specifying chromosome arms (DEFAULT: NULL)
#'
#' @return      A list, List, or *list, with pieces of x by chromosome arm
#'
#' @examples
#'
#' @export
#'
byChromArm <- function(x,
                       arms = NULL) { 

  if (is.null(arms) && 
      (!any(grepl("(GRCh(37|38)|hg(19|38))", unique(genome(x)))))) {
    stop("Your object's genome is odd (not hg19/GRCh37/38) and `arms` is NULL.")
  }

  if (is.null(arms)) {
    gome <- unique(genome(x))[1]
    arms <- paste0(gome, ".chromArm")
    data(list=arms, package="biscuiteer")
    arms <- get(arms)
  }

  ol <- findOverlaps(x, arms)
  rowRanges(x)$chromArm <- names(arms)[subjectHits(ol)]
  split(x, rowRanges(x)$chromArm) 

}
