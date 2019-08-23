#' wrapper fn to process large numbers of bigWig'ed CN segmentations
#' 
#' @param suffix    suffix of segmentation files (default is ".CN.hg19.bw")
#' @param filename  where to save the result (unsaved if NULL)
#' 
#' @return          a data.frame just like the one from grToSeg 
#' 
#' @export
bigWigsToSegs <- function(suffix=".CN.hg19.bw", filename=NULL) {
  bigWigs <- list.files(pattern=suffix)
  names(bigWigs) <- sub(suffix, "", bigWigs)
  grl <- GRangesList(lapply(bigWigs, import))
  grToSeg(grl, filename=filename)
}

