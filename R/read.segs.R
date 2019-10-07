#' Function to turn segmentations into GRanges with correction(s)
#'
#' If your segmentation file is not vaild, this probably won't work. A smarter
#' way to squash baselines would be with a mixture model. One bonus is that the
#' baselining is self-correcting - applying it to an already-corrected
#' segmentation will not do any harm and may even help. 
#'
#' @param file     The segmentation file to process
#' @param weight   Add weights for per-subject baselining? (DEFAULT: TRUE) 
#' @param by       If weighting, do so by marks or by width (DEFAULT: "marks")
#' @param from     Quantiles around the median to "squash" to baseline
#'                   (DEFAULT: c(.4, .6))
#' @param minNorm  Min log2ratio to consider as possible "normal" baseline
#'                   (DEFAULT: -0.5)
#' @param maxNorm  Max log2ratio to consider as possible "normal" baseline
#'                   (DEFAULT: +0.5)
#'
#' @return         A GRanges with corrections applied
#' 
#' @import GenomicRanges
#'
#' @export
#'
read.segs <- function(file,
                      weight = TRUE,
                      by = c("marks","width"), 
                      from = c(.4,.6),
                      minNorm = -0.5,
                      maxNorm = +0.5) {

  segs <- read.table(file, stringsAsFactors=FALSE, sep="\t")
  names(segs)[seq_len(4)] <- c("name","chrom","chromStart","chromEnd")
  names(segs)[ncol(segs) - 2] <- "width" # autocorrect below
  names(segs)[ncol(segs) - 1] <- "marks" # autocorrect below
  names(segs)[ncol(segs)] <- "score"
  if (weight) {
    by <- match.arg(by) 
    segs <- .addWeights(segs, by=by)
    segs <- .fixBaselines(segs, from=from, minNorm=minNorm, maxNorm=maxNorm)
    message("Baselines corrected.")
  }
  gr <- makeGRangesFromDataFrame(segs, keep.extra.columns=TRUE) 
  names(gr) <- gr$name
  return(gr)
}


# Helper function
.addWeights <- function(segs, by=c("marks","width")) {
  by <- match.arg(by)
  if (!by %in% names(segs)) {
    # in case the default does not / cannot work 
    segs$width <- segs$chromEnd - segs$chromStart
    by <- "width"
  }
  message("Weighting segments by ", by, "...", appendLF=FALSE)
  segs <- do.call(rbind, lapply(split(segs, segs$name), .addWeight, by=by))
  rownames(segs) <- segs$names
  message(" done.")
  return(segs)
}

# Helper function 
.addWeight <- function(seg, by=c("width","marks")) {
  by <- match.arg(by)  
  stopifnot(length(unique(seg$name)) == 1)
  if (by == "width") {
    seg$width <- seg$chromStart - seg$chromEnd  
    seg$weight <- seg$width / sum(seg$width)
  } else {
    seg$weight <- seg$marks / sum(seg$marks)
  }
  return(seg) 
}


# Helper function
.fixBaselines <- function(segs, from=c(.4, .6), minNorm=-0.5, maxNorm=+0.5) { 
  message("Finding baselines...", appendLF=FALSE)
  segs <- do.call(rbind, 
                  lapply(split(segs,segs$name), .fixBaseline, from=from,
                         minNorm=minNorm, maxNorm=maxNorm))
  rownames(segs) <- segs$names
  message(" done.")
  return(segs)
}


# Helper function 
#' @importFrom stats quantile
.fixBaseline <- function(seg, from=c(.4, .6), minNorm=-0.5, maxNorm=+0.5) { 
  baseline <- .findBaseline(seg, minNorm=minNorm, maxNorm=maxNorm)
  higher <- quantile(seg$score, prob=max(from), weight=seg$weight, na.rm=TRUE)
  lower <- quantile(seg$score, prob=min(from), weight=seg$weight, na.rm=TRUE)
  probablyUnchanged <- which(seg$score >= lower & seg$score <= higher)
  seg[probablyUnchanged, "score"] <- baseline
  seg$score <- seg$score - baseline
  return(seg)
}


# Helper function
#' @importFrom stats median
.findBaseline <- function(seg, minNorm=-0.5, maxNorm=0.5) { 
  nearNormal <- subset(seg, score >= minNorm & score <= maxNorm)
  median(nearNormal$score, weight=nearNormal$weight, na.rm=TRUE)
}


# Helper function
# TODO: Work in progress
.assignUncertainty <- function(seg, plot=FALSE) { 
  # correct for baseline (median correction; only works for mostly-NK tumors)
  seg$score <- seg$score - .findBaseline(seg)

  # generate Rle representations of score X marks, post-shift  

  # unlist into an "unweighted" vector of weighted scores

  # run mclust on that ?

  # assign segments in the "baseline" cluster a log-ratio of 0. 
  
  # plot uncertainty, if requested 

  # return the corrected segments. 
}


# Helper function
.write.segs <- function(gr, file=NULL) { 
  names(gr) <- NULL 
  segs <- as.data.frame(gr)
  segs <- segs[, c("name","seqnames","start","end","score")] 
  if (!is.null(file)) {
    write.table(segs, sep="\t", quote=FALSE, row.names=FALSE,
                col.names=FALSE, file=file) 
  }
  return(segs) 
}
