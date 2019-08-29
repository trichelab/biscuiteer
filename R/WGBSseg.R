#' Copy number segmentation of GC- & mapping-corrected WGBS coverage
#'
#' FIXME: use some form of smoothing to glob segments together and THEN filter 
#'        on length (as opposed to right now where it is done rather roughly)
#' FIXME: use findOverlaps(resize(segs, width(segs) + minWidth, fix="center"))
#' 
#' @param cn    	binned, corrected read depths from correctBsSeqCoverage()
#' @param minAbs	minimum absolute value of a segment to keep it (0.1)
#' @param minWidth      minimum width of segment (NOT A GAP) to keep (25000)
#' @param ...           other parameters passed to fastseg
#'
#' @return              a GRanges (hg19 seqinfo)
#'
#' @import fastseg
#' @import Homo.sapiens
#' 
#' @examples
#' 
#'   library(fastseg)
#'   library(HMMcopy) 
#'   library(rtracklayer)
#'   files <- list.files(patt="P01-022-.*covg.hg19.bw") 
#'   names(files) <- sapply(strsplit(files, "(\\.|\\-)"), `[`, 3)
#'   covgs <- sapply(files, import)
#'
#'   RMS1 <- correctBsSeqCoverage(covgs$tumor1, covgs$normal)
#'   RMS1seg <- WGBSseg(RMS1)
#'   plotSegs(RMS1, RMS1seg, "chr10")
#'
#'   RMS2 <- correctBsSeqCoverage(covgs$tumor2, covgs$normal)
#'   RMS2seg <- WGBSseg(RMS2) 
#'   plotSegs(RMS2, RMS2seg, "chr10")
#' 
#' @export
WGBSseg <- function(cn, minAbs=0.3, minWidth=20000, ...) {

  if (is.null(attr(cn, "corrected"))) {
    stop("You need to run correctBsSeqCoverage(tumor, normal) first...")
  }
  mcols(cn) <- mcols(cn)[,"copy"]
  message("Segmenting...")
  res <- suppressMessages(fastseg(cn, ...)) 
  message("Verifying seqinfo()...")
  if (anyNA(seqlengths(res))) {
    # needed by gaps() in the minAbs filtering code below 
    seqinfo(res) <- seqinfo(Homo.sapiens)[seqlevels(res)]
  }
  res$state <- 3
  res$median <- res$seg.mean # gross...

  mergeGaps <- function(res, gap, minWidth=0) {
    if (!any(width(gap) < minWidth)) return(gap[0])
    gap <- subset(gap, width(gap) < minWidth)
    gap <- resize(gap, width(gap) + minWidth, fix="center")
    gap <- subsetByOverlaps(gap, res)
    gapResOlaps <- findOverlaps(gap, res)
    gap$state <- 3
    gap$median <- tapply(res$median[subjectHits(gapResOlaps)],
                         queryHits(gapResOlaps), median)
    gap <- unique(subset(gap, abs(gap$median) >= minAbs))
    gap <- resize(gap, width(gap)-minWidth, fix="center")
    gap$state[which(gap$median >= minAbs)] <- 4
    gap$state[which(gap$median >= 0.6)] <- 5
    gap$state[which(gap$median >= 1.2)] <- 6
    gap$state[which(gap$median <= -1 * minAbs)] <- 2
    gap$state[which(gap$median <= -1.2)] <- 1
    mcols(gap) <- mcols(gap)[, c("state","median")]
    return(gap)
  }

  message("Merging gaps less than ", minWidth, " bases wide...")
  if (any(abs(res$seg.mean) >= minAbs & width(res) >= minWidth)) { 

    # effetively purge most same-sign gaps less than minWidth across:
    res <- subset(res, width(res) >= minWidth)
    res$state[which(res$seg.mean >= minAbs)] <- 4
    res$state[which(res$seg.mean >= 0.6)] <- 5
    res$state[which(res$seg.mean >= 1.2)] <- 6
    res$state[which(res$seg.mean <= -1 * minAbs)] <- 2
    res$state[which(res$seg.mean <= -1.2)] <- 1
    mcols(res) <- mcols(res)[, c("state","median")]
    res <- subset(res, state != 3)
    gap <- subset(gaps(res), strand=="*")
    gap$state <- 3
    gap$median <- 0
    res <- c(res, mergeGaps(res, gap, minWidth))
    gap <- subset(gaps(res), strand=="*")
    gap$state <- 3
    gap$median <- 0
    res <- sort(c(gap, res))
  }
  attr(res, "param") <- list(...)
  res$score <- res$median
  return(res)

}

