#' Process many bigWig'ed CN segmentations
#'
#' Wrapper function to grToSeg to make processing large numbers of bigWigs easier
#'
#' @param bwList    bigWig files to read in
#' @param suffix    Suffix of segmentation files (DEFAULT: ".CN.hg19.bw")
#' @param filename  Where to save the result - unsaved if NULL (DEFAULT: NULL)
#'
#' @return          A data.frame like the one returned from grToSeg
#'
#' @examples
#'
#' @export
#'
bigWigsToSegs <- function(bwList,
                          suffix = ".CN.hg19.bw",
                          filename = NULL) {
  bigWigs <- bwList
  names(bigWigs) <- sub(suffix, "", bigWigs)
  grl <- GRangesList(lapply(bigWigs, import))
  grToSeg(grl, filename=filename)
}

