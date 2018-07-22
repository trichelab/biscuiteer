#' an idiotically simple parallelization step that I've been using for years
#'
#' Since this is an idiotically simple function, it uses the idiot's favorite 
#' genome assembly, which is hg19*. If you're going to be a smart person, you
#' will need to either supply another argument `arms` with those ranges, or 
#' send in a patch so we can properly support arbitrary genomes (or maybe push
#' it up the line to GenomeInfoDb for this super generic feature, hint hint).
#'
#' This function just splits an object by chromosome arm, which tends to make
#' parallelization MUCH easier, as cross-arm dependencies are unusual, and the 
#' larger chromosomes can thus be split across processes or machines without 
#' worrying much about data starvation for processes on smaller chromosomes.
#' 
#' * Note that it's OK to be smart and use hg19_rCRSchrm,
#'   because mitochondrial chromosomes don't have arms. 
#'
#' @param x     anything with a GRanges in it: BSseq, SummarizedExperiment, etc.
#' @param arms  (optional) another GRanges, but just specifying chromosome arms.
#'
#' @return      a list, List, or *List, with pieces of x by chromosome arm.
#' 
#' @export
byChromArm <- function(x, arms=NULL) { 

  if (is.null(arms) && (!any(grepl("(GRCh37|hg19)", unique(genome(x)))))) {
    stop("Your object's genome isn't labeled as GRCh37/hg19 and `arms` is NULL")
  }

  stop("Not quite finished") 

}
