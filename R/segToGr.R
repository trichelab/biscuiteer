#' Import a segmentation file into GRanges object
#'
#' Reverse of grToSeg
#'
#' @param seg     The .seg filename
#' @param genome  Genome against which segments were annotated (DEFAULT: "hg19")
#' @param name    .seg file column to use as $name metadata (DEFAULT: "ID")
#' @param score   .seg file column to use as $score metadata
#'                  (DEFAULT: "seg.mean")
#' 
#' @return        A GRanges object
#'
#' @import GenomicRanges
#'
#' @seealso grToSeg
#'
#' @examples
#'
#' @export
#'
segToGr <- function(seg,
                    genome = "hg19",
                    name = "ID",
                    score = "seg.mean") { 
  segs <- read.table(seg, header=TRUE, sep="\t")
  names(segs) <- base::sub(name, "name", names(segs))
  names(segs) <- base::sub(score, "score", names(segs))
  firstTwo <- c("score", "name")
  segs <- segs[, c(firstTwo, setdiff(names(segs), firstTwo))]
  if (nrow(segs) > 0) {
    gr <- sort(makeGRangesFromDataFrame(segs, keep=TRUE))
    names(gr) <- gr$name
    return(gr)
  } else { 
    return(GRanges())
  }
}

