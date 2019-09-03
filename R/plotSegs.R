#' Plot WGBScopy results
#'
#' Modified from HMMcopy
#'
#' @param cn    Binned, corrected read depths from correctBsSeqCoverage
#' @param segs  Segmentation result from WGBScopy
#' @param chr   Chromosome to display (DEFAULT: "chr22")
#' @param name  Name for the plot (DEFAULT: "Segmented")
#' @param ...   Other arguments to pass on to plot
#'
#' @import GenomicRanges
#' @import HMMcopy
#'
#' @examples
#'
#' @export
#'
plotSegs <- function(cn,
                     segs = NULL,
                     chr = "chr22",
                     name = "Segmented",
                     ...) {

  # accept a list from TN_WGBS_CNA()
  if (is.null(segs) & is.list(cn)) {
    segs <- cn$segs
    cn <- cn$cn
  }
  cn$state <- 3 
  if (!chr %in% seqnames(cn)) stop(paste(chr, "is not among seqnames(cn)"))
  olaps <- findOverlaps(cn, segs)
  cn$state[queryHits(olaps)] <- segs$state[subjectHits(olaps)]
  cols <- c("darkblue","blue","gray80","#FFCCCC","red","darkred")
  range <- quantile(cn$copy, na.rm=TRUE, prob=c(0.00001, 0.99999))
  a <- subset(cn, seqnames(cn) == chr)
  b <- subset(segs, seqnames(segs) == chr)
  plot(start(a), a$copy, col=cols[as.numeric(as.character(a$state))], 
       ylim=range, ylab="Tumor copy number", xlab="Chromosome position", 
       main=paste(name, chr), ...)
  for (k in seq_along(b)) {
    lines(c(start(b)[k], end(b)[k]), rep(b$median[k], 2), lwd=3, col="green")
  } 
  legend("topleft", c("HomDel", "Del", "Neutral", "Gain", "Amp", "BigAmp"),
         fill=cols, horiz=TRUE, bty="n", cex=1)
}

