#' Guess ages using Horvath-style 'clock' models
#'
#' See Horvath, Genome Biology, 2013 for more information
#'
#' Note: the accuracy of the prediction will increase or decrease depending on
#' how various hyper-parameters are set by the user. This is NOT a hands-off 
#' procedure, and the defaults are only a starting point for exploration. It
#' will not be uncommon to tune `padding`, `minCovg`, and `minSamp` for each
#' WGBS or RRBS experiment (and the latter may be impacted by whether dupes are
#' removed prior to importing data). Consider yourself forewarned. In the near
#' future we may add support for arbitrary region-coefficient inputs and result
#' transformation functions, which of course will just make the problems worse.
#' 
#' Also, please cite the appropriate papers for the Epigenetic Clock(s) you use:
#'
#' For the 'horvath' or 'horvathshrunk' clocks, cite Horvath,
#' Genome Biology 2013.
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
#' @param bsseq    A bsseq object (must have assays named `M` and `Cov`)
#' @param model    Which model ("horvath", "horvathshrunk", "hannum",
#'                   "skinandblood")
#' @param padding  How many bases +/- to pad the target CpG by (DEFAULT: 15)
#' @param useENSR  Use ENSEMBL regulatory region bounds instead of CpGs
#'                   (DEFAULT: FALSE)
#' @param useHMMI  Use HMM CpG island boundaries instead of padded CpGs
#'                   (DEFAULT: FALSE)
#' @param minCovg  Minimum regional read coverage desired to estimate 5mC
#'                   (DEFAULT: 5)
#' @param impute   Use k-NN imputation to fill in low-coverage regions?
#'                   (DEFAULT: FALSE)
#' @param minSamp  Minimum number of non-NA samples to perform imputation
#'                   (DEFAULT: 5)
#' @param genome   Genome to use as reference, if no genome(bsseq) is set
#'                   (DEFAULT: NULL)
#' @param dropBad  Drop rows/cols with > half missing pre-imputation?
#'                   (DEFAULT: FALSE)
#' @param ...      Arguments to be passed to impute.knn, such as rng.seed
#'  
#' @return         A list: call, methylation estimates, coefs, age (estimates)
#'
#' @import impute
#' @importFrom methods as
#'
#' @examples
#' 
#'   shuf_bed <- system.file("extdata", "MCF7_Cunha_chr11p15_shuffled.bed.gz",
#'                           package="biscuiteer")
#'   orig_bed <- system.file("extdata", "MCF7_Cunha_chr11p15.bed.gz",
#'                           package="biscuiteer")
#'   shuf_vcf <- system.file("extdata",
#'                           "MCF7_Cunha_shuffled_header_only.vcf.gz",
#'                           package="biscuiteer")
#'   orig_vcf <- system.file("extdata",
#'                           "MCF7_Cunha_header_only.vcf.gz",
#'                           package="biscuiteer")
#'   bisc1 <- read.biscuit(BEDfile = shuf_bed, VCFfile = shuf_vcf,
#'                         merged = FALSE)
#'   bisc2 <- read.biscuit(BEDfile = orig_bed, VCFfile = orig_vcf,
#'                         merged = FALSE)
#'
#'   comb <- unionize(bisc1, bisc2)
#'   ages <- WGBSage(comb, "horvath")
#'
#' @export
#'
WGBSage <- function(bsseq,
                    model = c("horvath","horvathshrunk",
                              "hannum","skinandblood"),
                    padding = 15,
                    useENSR = FALSE,
                    useHMMI = FALSE, 
                    minCovg = 5,
                    impute = FALSE,
                    minSamp = 5,
                    genome = NULL, 
                    dropBad = FALSE,
                    ...) { 

  # sort out assemblies
  g <- unique(genome(bsseq))
  if (is.null(g)) g <- genome 
  if (!g %in% c("hg19","GRCh37","hg38","GRCh38")) {
    message("genome(bsseq) is set to ",unique(genome(bsseq)),
            ", which is unsupported.")
    stop("Provide a `genome` argument, ",
         "or set genome(bsseq) manually, to proceed.")
  }

  # get the requested model (a List with regions, intercept, and a cleanup fn
  model <- match.arg(model)
  clock <- getClock(model=model, padding=padding, genome=g,
                    useENSR=useENSR, useHMMI=useHMMI)
 
  # assess coverage (since this affects the precision of estimates)
  message("Assessing coverage across age-associated regions...") 
  covgWGBSage <- getCoverage(bsseq, regions=clock$gr, what="perRegionTotal")
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
  methWGBSage <- getMeth(bsseq, regions=clock$gr, type="raw", what="perRegion")
  rownames(methWGBSage) <- as.character(clock$gr)
  methWGBSage[covgWGBSage < minCovg] <- NA
  methWGBSage <- as(methWGBSage, "matrix") # 353 x ncol(bsseq) is not too huge

  # drop samples if needed
  droppedSamples <- c()
  droppedRegions <- c()
  if (dropBad) {
    thresh <- ceiling(length(clock$gr) / 2)
    droppedSamples <- names(which(colSums(is.na(methWGBSage)) >= thresh))
    keptSamples <- setdiff(colnames(bsseq), droppedSamples)
    thresh2 <- ceiling(ncol(bsseq) / 2)
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
