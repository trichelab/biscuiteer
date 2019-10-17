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
#' @importFrom utils data
#' @import SummarizedExperiment
#'
#' @examples
#'
#'   orig_bed <- system.file("extdata", "MCF7_Cunha_chr11p15.bed.gz",
#'                           package="biscuiteer")
#'   orig_vcf <- system.file("extdata", "MCF7_Cunha_header_only.vcf.gz",
#'                           package="biscuiteer")
#'   bisc <- readBiscuit(BEDfile = orig_bed, VCFfile = orig_vcf,
#'                       merged = FALSE)
#'
#'   reg <- GRanges(seqnames = rep("chr11",5),
#'                  strand = rep("*",5),
#'                  ranges = IRanges(start = c(0,2.8e6,1.17e7,1.38e7,1.69e7),
#'                                   end= c(2.8e6,1.17e7,1.38e7,1.69e7,2.2e7))
#'                  )
#'   names(reg) <- as.character(reg)
#'
#'   arms <- byChromArm(bisc, arms = reg)
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
