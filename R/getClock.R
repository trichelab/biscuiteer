#' Retrieve 'epigenetic clock' models
#'
#' Biscuiteer supports several 'epigenetic clock' models. This function
#' retrieves the various models.
#'
#' The remapped coordinates for the Horvath (2012) and Hannum (2013) clocks, 
#' along with shrunken Horvath (2012) and improved Horvath (2018) models, are
#' provided as part of biscuiteer (visit inst/extdata/clocks.R to find out how)
#' along with some functionality to make them more usable in RRBS/WGBS data of 
#' varying coverage along varying genomes. For example, the HMM-based CpG island
#' model introduced by Wu (2010) can be used to assign to within-island features
#' the methylation rate of their associated island, and ENSEMBL regulatory build
#' features (ENSR features, for short) such as CTCF binding sites can have their
#' coordinates substituted for the default padded boundaries of a feature. 
#' 
#' The net result of this process is that, while the default settings simply 
#' swap in a 30-bp stretch centered on the selected clock's CpG (and/or CpH) 
#' loci, add the intercept, and ship out the model, much more flexibility is 
#' available to the user. This function provides a single point for tuning of
#' such options in the event that defaults don't work well for a user. 
#' 
#' The precedence of options is as follows: 
#' \enumerate{
#'   \item If a feature has neither ENSR nor HMMI IDs, it is padded
#'         (only) +/- bp.
#'   \item If it has an HMMI but not ENSR ID or ENSR==FALSE, the HMM island
#'         is used.
#'   \item If a feature has an ENSR ID, and ENSR==TRUE, the ENSR feature
#'         is used.
#' }
#'
#' If a feature has both an ENSR ID and an HMMI ID, and both options are TRUE, 
#' then the ENSR start and end coordinates will take precedence over its HMMI. 
#'
#' The above shenanigans produce the GRanges object returned as `gr` in a List.
#' The `intercept` value returned with the model is its fixed (B0) coefficient.
#' The `cleanup` function returned with the model transforms its raw output. 
#' 
#' @param model    One of "horvath", "horvathshrunk", "hannum", or
#'                   "skinandblood"
#' @param padding  How many base pairs (+/-) to expand a feature's footprint
#'                   (DEFAULT: 15)
#' @param genome   One of "hg19", "GRCh37", "hg38", or "GRCh38"
#'                   (DEFAULT: "hg19")
#' @param useENSR  Substitute ENSEMBL regulatory feature boundaries?
#'                   (DEFAULT: FALSE) 
#' @param useHMMI  Substitute HMM-based CpG island boundaries?
#'                   (DEFAULT: FALSE) 
#'
#' @return         a List with elements `model`, `gr`, `intercept`,
#'                   and `cleanup`
#' 
#' @examples
#'
#' clock <- getClock(model="horvathshrunk", genome="hg38")
#'
#' @export
#'
getClock <- function(model = c("horvath","horvathshrunk",
                               "hannum","skinandblood"),
                     padding = 15,
                     genome = c("hg19","hg38","GRCh37","GRCh38"),
                     useENSR = FALSE,
                     useHMMI = FALSE) {

  model <- match.arg(model) 
  genome <- match.arg(genome)
  g <- sub("hg37","hg19", sub("GRCh", "hg", genome))
  grcols <- paste0(g, c("chrom","start","end","HMMI","ENSR"))
  data(clocks, package="biscuiteer")
    
  clock <- subset(clocks[, c(model, grcols)], !is.na(clocks[, model]))
  gr <- sort(makeGRangesFromDataFrame(clock[-1,], 
                                      seqnames.field=grcols[1], 
                                      start.field=grcols[2], 
                                      end.field=grcols[3], 
                                      keep.extra.columns=TRUE))
  names(mcols(gr)) <- c("score", "HMMI", "ENSR")

  # expand the ranges by adding padding bases around the target loci? 
  if (padding > 0) gr <- trim(resize(gr, 2 * padding, fix="center"))
  
  # use islands?
  if (useHMMI) {
    hmmdat <- paste0("HMM_CpG_islands.", g)
    data(list=hmmdat)
    HMMIs <- get(hmmdat)
    HMMIed <- names(subset(gr, !is.na(HMMI)))
    ranges(gr)[HMMIed] <- ranges(HMMIs[gr[HMMIed]$HMMI])
    gr$HMMI <- NULL # not really necessary post-expansion
  } 

  # use ENSRs?
  if (useENSR) {
    ensrdat <- paste0("ENSR_subset.", g)
    data(list=ensrdat)
    ENSRs <- get(ensrdat)
    ENSRed <- names(subset(gr, !is.na(ENSR)))
    ranges(gr)[ENSRed] <- ranges(ENSRs[gr[ENSRed]$ENSR])
  }

  # fix seqinfo
  genome(gr) <- g
  seqdat <- paste0("seqinfo.", g)
  data(list=seqdat)
  seqinfo(gr) <- get(seqdat)[seqlevels(gr)] 

  # switch naming styles if GRCh37 or GRCh38
  genome(gr) <- genome
  if (genome %in% c("GRCh37","GRCh38")) seqlevelsStyle(gr) <- "NCBI" 

  # swap in the right cleanup fn (all but Hannum use an exponential transform)
  cleanup <- switch(model, hannum=base::identity, biscuiteer::fixAge)

  # done
  res <- List(model=model, gr=gr, intercept=clock[1,1], cleanup=cleanup)
  return(res)

}
