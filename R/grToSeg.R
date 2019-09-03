#' Dump GRanges to segmented data data.frame
#'
#' Output data.frame can be written to a .seg file if supplied with filename
#' input argument
#'
#' @param gr        A GRanges or GRangesList to dump to .seg file
#' @param filename  Where to save the result - unsaved if NULL (DEFAULT: NULL)
#' @param minAbs    Minimum absolute gain/loss cutoff (DEFAULT: NULL)
#'
#' @return          A data.frame with columns:
#'                    (ID, chrom, loc.start, loc.end, num.mark, seg.mean)
#'
#' @import GenomicRanges
#'
#' @seealso segToGr
#'
#' @examples
#'
#' @export
#'
grToSeg <- function(gr,
                    filename = NULL,
                    minAbs = NULL) { 
  if (is(gr, "GRangesList")) {
    seglist <- lapply(gr, grToSeg, filename=NULL)
    for (i in names(seglist)) seglist[[i]][, "ID"] <- i
    segs <- do.call(rbind, seglist)
  } else { 
    if (!is.null(minAbs)) gr <- subset(gr, abs(gr$score) >= minAbs)
    if (is.null(names(gr))) names(gr) <- paste0("segment", seq_along(gr))
    segs <- cbind("ID"=names(gr), as.data.frame(gr))
    names(segs) <- sub("^seqnames$", "chrom", names(segs))
    names(segs) <- sub("^start$", "loc.start", names(segs))
    names(segs) <- sub("^end$", "loc.end", names(segs))
    names(segs) <- sub("^score$", "seg.mean", names(segs))
    segs$num.mark <- round(segs$width / 1000) # 1kb bins per WGBS CNA
    segs <- segs[, c("ID","chrom","loc.start","loc.end","num.mark","seg.mean")]
  }
  rownames(segs) <- NULL
  if (!is.null(filename)) {
    write.table(segs, file=filename, quote=FALSE, sep="\t", row.names=FALSE)
  }
  return(segs)
}

