#' tweaked version of HMMcopy's GC bias and mappability function
#'
#' @param x           binned hg19 coverage with GC content and mappability cols
#' @param mappability minimum mappability to use bins as examplars for others
#' @param samplesize  lowess sample size (points to sample; default 50000)
#' @param verbose     be chatty? (TRUE)
#'
#' @return            corrected tumor and normal read counts
#'
#' @import GenomicRanges
#' @import HMMcopy
#' 
#' @export
correctDepth <- function(x, mappability=0.9, samplesize=50000, verbose=TRUE) {

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
