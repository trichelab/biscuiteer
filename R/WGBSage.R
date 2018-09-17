#' Guess ages using various Horvath-style (Genome Biology, 2012) 'clock' models
#'
#' Note: the accuracy of the prediction will increase or decrease depending on
#' how various [hyper]parameters are set by the user. This is NOT a hands-off 
#' procedure, and the defaults are only a starting point for exploration. It
#' will not be uncommon to tune `padding`, `minCovg`, and `minSamp` for each
#' WGBS or RRBS experiment (and the latter may be impacted by whether dupes are
#' removed prior to importing data). Consider yourself forewarned. In the near
#' future we may add support for arbitrary region-coefficient inputs and result
#' transformation functions, which of course will just make the problems worse.
#' 
#' Also, please cite the appropriate papers for the Epigenetic Clock(s) you use:
#'
#' For the 'horvath' or 'shrunk' clocks, cite Horvath, Genome Biology 2012.
#' For the 'hannum' clock, cite Hannum et al, Molecular Cell 2013. 
#' For the 'skinandblood' clock, cite Horvath et al, Aging 2018. 
#' 
#' Last but not least, keep track of the parameters YOU used for YOUR estimates.
#' The `call` element in the returned list of results is for this exact purpose.
#' If you need to recover the GRanges object used to average (or impute) DNAme 
#' values for the model, try `as.character(rownames(result$meth))` on a result.
#' The coefficients for each of these regions are stored in result$coefs, and 
#' the age estimates are stored in result$age (named, in case dropBad == TRUE). 
#' 
#' @param   x       a BSseq object (must have assays named `M` and `Cov`)
#' @param   model   which model ("horvath","shrunk","hannum","skinandblood")
#' @param   padding how many bases +/- to pad the target CpG by (default is 15)
#' @param   useENSR use ENSEMBL regulatory region bounds instead of CpGs (FALSE)
#' @param   useHMMI use HMM CpG island boundaries instead of padded CpGs (FALSE)
#' @param   minCovg minimum regional read coverage desired to estimate 5mC (5)
#' @param   impute  use k-NN imputation to fill in low-coverage regions? (FALSE)
#' @param   minSamp minimum number of non-NA samples to perform imputation (5)
#' @param   genome  genome to use as reference, if no genome(x) is set (NULL) 
#' @param   dropBad drop rows/cols with > half missing pre-imputation? (FALSE) 
#' @param   ...     arguments to be passed to impute.knn, such as rng.seed
#'  
#' @return          a list: call, methylation estimates, coefs, age (estimates)
#'
#' @import  impute
#' 
#' @export
WGBSage <- function(x, model=c("horvath","shrunk","hannum","skinandblood"),
                    padding=15, useENSR=FALSE, useHMMI=FALSE, 
                    minCovg=5, impute=FALSE, minSamp=5, genome=NULL, 
                    dropBad=FALSE, ...) { 

  # sort out assemblies
  g <- unique(genome(x))
  if (is.null(g)) g <- genome 
  if (!g %in% c("hg19","GRCh37","hg38","GRCh38")) {
    message("genome(x) is set to ",unique(genome(x)),", which is unsupported.")
    stop("Provide a `genome` argument, or set genome(x) manually, to proceed.")
  }

  # get the requested model (a List with regions, intercept, and a cleanup fn
  model <- match.arg(model)
  clock <- getClock(model=model, padding=padding, genome=g,
                    useENSR=useENSR, useHMMI=useHMMI)
 
  # assess coverage (since this affects the precision of estimates)
  message("Assessing coverage across age-associated regions...") 
  covgWGBSage <- getCoverage(x, regions=clock$gr, what="perRegionTotal")
  rownames(covgWGBSage) <- names(clock$gr)
  NAs <- which(as.matrix(covgWGBSage) < minCovg, arr.ind=TRUE)

  # for sample/region dropping
  subM <- rep(FALSE, nrow(covgWGBSage))
  names(subM) <- rownames(covgWGBSage)

  # tabulate for above
  if (nrow(NAs) > 0) {
    # either way, probably a good idea to fix stuff 
    message("You have NAs. Change `padding` (",padding,"), `minCovg` (",
            minCovg,"), `useHMMI`, and/or `useENSR`.")
  } else { 
    message("All regions in all samples appear to be sufficiently covered.") 
  }

  # extract regional methylation estimates (and set < minCovg sites to NA) 
  methWGBSage <- getMeth(x, regions=clock$gr, type="raw", what="perRegion")
  rownames(methWGBSage) <- as.character(clock$gr)
  methWGBSage[covgWGBSage < minCovg] <- NA
  methWGBSage <- as(methWGBSage, "matrix") # 353 x ncol(x) is not too huge

  # drop samples if needed
  droppedSamples <- c()
  droppedRegions <- c()
  if (dropBad) {
    thresh <- ceiling(length(clock$gr) / 2)
    droppedSamples <- names(which(colSums(is.na(methWGBSage)) >= thresh))
    keptSamples <- setdiff(colnames(x), droppedSamples)
    thresh2 <- ceiling(ncol(x) / 2)
    droppedRegions <- names(which(rowSums(is.na(methWGBSage)) >= thresh2))
    keptRegions <- setdiff(names(clock$gr), droppedRegions)
    methWGBSage <- methWGBSage[keptRegions, keptSamples] 
  }
  
  # impute, if requested, any sites with insufficient coverage
  if (impute) methWGBSage <- impute.knn(methWGBSage, k=minSamp, ...)$data
  keep <- (rowSums2(is.na(methWGBSage)) < 1)
  if (!all(keep)) methWGBSage <- methWGBSage[which(keep), ] 

  names(clock$gr) <- as.character(granges(clock$gr))
  coefs <- clock$gr[rownames(methWGBSage)]$score
  names(coefs) <- rownames(methWGBSage)
  agePredRaw <- (clock$intercept + (t(methWGBSage) %*% coefs))
  res <- list(call=sys.call(), 
              droppedSamples=droppedSamples,
              droppedRegions=droppedRegions,
              intercept=clock$intercept,
              meth=methWGBSage, 
              coefs=coefs,
              age=clock$cleanup(agePredRaw))
  return(res)
    
}
