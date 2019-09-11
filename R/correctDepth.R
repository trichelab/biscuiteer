#' Tweaked version of HMMcopy's GC bias and mappability function
#'
#' Also used by correctBsSeqCoverage
#'
#' @param x            Binned hg19 coverage with GC content and mappability
#'                       columns
#' @param mappability  Mappability threshold below which points are ignored
#'                       while creating the correction curve (value between 0
#'                       and 1)
#' @param samplesize   Number of points sampled during LOESS fitting
#'                       (DEFAULT: 50 000)
#' @param verbose      Print extra statements? (DEFAULT: TRUE)
#'
#' @return             Corrected tumor and normal read counts
#'
#' @import GenomicRanges
#' @import HMMcopy
#'
#' @seealso correctBsSeqCoverage
#' @seealso HMMcopy::correctReadcount
#'
#' @examples
#'
#' @export
#'
correctDepth <- function(x,
                         mappability = 0.9,
                         samplesize = 50000,
                         verbose = TRUE) {

  if (!all( c("reads","gc","map") %in% names(mcols(x)))) {
    stop("Missing one of required columns: reads, gc, map")
  }

  if(verbose) message("Applying filter on data...") 
  x$valid <- TRUE
  x$valid[x$reads <= 0 | x$gc < 0] <- FALSE
  x$ideal <- TRUE
  routlier <- 0.01
  range <- quantile(x$reads[x$valid], prob=c(0, 1-routlier), na.rm=TRUE)
  doutlier <- 0.001
  domain <- quantile(x$gc[x$valid], prob=c(doutlier, 1-doutlier), na.rm=TRUE)
  x$ideal[!x$valid | x$map < mappability | x$reads <= range[1] |
    x$reads > range[2] | x$gc < domain[1] | x$gc > domain[2]] <- FALSE

  if (verbose) message("Correcting for GC bias...")
  set <- which(x$ideal)
  select <- sample(set, min(length(set), samplesize))
  rough <- loess(x$reads[select] ~ x$gc[select], span = 0.03)
  i <- seq(0, 1, by = 0.001)
  final = loess(predict(rough, i) ~ i, span = 0.3)
  x$cor.gc <- x$reads / predict(final, x$gc)

  if (verbose) message("Correcting for mappability bias...")
  coutlier <- 0.01
  range <- quantile(x$cor.gc[which(x$valid)], prob=c(0, 1-coutlier), na.rm=TRUE)
  set <- which(x$cor.gc < range[2])
  select <- sample(set, min(length(set), samplesize))
  final <- approxfun(lowess(x$map[select], x$cor.gc[select]))
  x$cor.map <- x$cor.gc / final(x$map)
  x$copy <- x$cor.map
  x$copy[x$copy <= 0] = NA
  x$copy <- log(x$copy, 2)
  x <- keepSeqlevels(x, unique(seqnames(x)))
  attr(x, "corrected") <- TRUE 
  return(x)
}
