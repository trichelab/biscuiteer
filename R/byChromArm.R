#' an idiotically simple parallelization step that I've been using for years
#'
#' This function just splits an object by chromosome arm, which tends to make
#' parallelization MUCH easier, as cross-arm dependencies are unusual, and the 
#' larger chromosomes can thus be split across processes or machines without 
#' worrying much about data starvation for processes on smaller chromosomes.
#' 
#' @param x     anything with a GRanges in it: BSseq, SummarizedExperiment, etc.
#' @param arms  (optional) another GRanges, but just specifying chromosome arms.
#'
#' @return      a list, List, or *List, with pieces of x by chromosome arm.
#' 
#' @export
byChromArm <- function(x, arms=NULL) { 

  if (is.null(arms) && 
      (!any(grepl("(GRCh(37|38)|hg(19|38))", unique(genome(x)))))) {
    stop("Your object's genome is odd (not hg19/GRCh37/38) and `arms` is NULL.")
  }

  if (is.null(arms)) {
    gome <- unique(genome(x))[1]
    arms <- paste0(gome, ".chromArm")
    data(list=arms, package="biscuitEater")
    arms <- get(arms)
  }

  ol <- findOverlaps(x, arms)
  rowRanges(x)$chromArm <- names(arms)[subjectHits(ol)]
  split(x, rowRanges(x)$chromArm) 

}
