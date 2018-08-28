#' bin read coverage to simplify and improve CNA calling 
#' 
#' helper function for binning .covg.hg19.bw bigWigs (see below for note)
#' nb. it is best to extract these from Biscuit BAMs via genomeCoverageBed
#' FIXME: figure out how to estimate the most likely GCbias ~ DNAme linkage
#'
#' @param covg      the coverage track (a GRanges) or a bigWig filename 
#' @param bins      the bins to summarize over (use tileGenome to generate)
#' @param tilewidth if is.null(bins), what width should be used? (1000bp)
#' @param which     if covg is a filename, what regions should be imported?
#'
#' @return          binned read counts
#'
#' @import GenomicRanges
#' @import Homo.sapiens
#' @import rtracklayer
#' 
#' @export
binCoverage <- function(covg, bins=NULL, tilewidth=1000, which=NULL) {

  if (is(covg, "character")) covg <- import(covg)
  if (is.null(bins)) {
    bins <- tileGenome(seqlengths(Homo.sapiens), 
                       tilewidth=tilewidth, 
                       cut.last=TRUE)
  }
  bins$score <- 0 
  olaps <- findOverlaps(covg, bins, type="within")
  summed <- tapply(covg$score[queryHits(olaps)], subjectHits(olaps), sum)
  bins$score[as.integer(names(summed))] <- summed 
  attr(bins, "binned") <- TRUE

  if (!is.null(which)) bins <- subsetByOverlaps(bins, which)
	
  return(bins)

}
