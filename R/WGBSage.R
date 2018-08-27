#' Guess specimens' ages by a modified Horvath (Genome Biology, 2013) algorithm
#'
#' Note: the accuracy of the prediction will increase or decrease depending on
#' how the various parameters are set by the user. This is NOT a hands-off 
#' procedure, and the defaults are only a starting point for exploration. It
#' will not be uncommon to tune `pad`, `minCovg`, and `minSamp` for each WGBS
#' or [e]RRBS experiment (and the latter may be impacted by whether dupes are 
#' removed prior to importing data). Please consider yourself forewarned!
#' 
#' @param   x       a BSseq object (must have assays named `M` and `Cov`)
#' @param   pad     how many bases to pad the target CpG by (default is 15)
#' @param   minCovg minimum regional read coverage desired to estimate 5mC% (5)
#' @param   impute  use k-NN imputation to fill in low-coverage regions? (TRUE) 
#' @param   minSamp minimum number of non-NA samples to perform imputation (5)
#' @param   shrink  use shrunken version of Horvath's coefs & intercept (FALSE)
#' @param   useENSR use ENSEMBL regulatory region bounds instead of CpGs (FALSE)
#' @param   genome  genome to use as reference, if no genome(x) is set (NULL) 
#' @param   ...     arguments passed to impute.knn, such as rng.seed
#'  
#' @return          a list: age estimates, meth estimates, and parameters
#'
#' @import  impute
#' 
#' @export
WGBSage <- function(x, pad=15, minCovg=5, impute=TRUE, minSamp=5, shrink=FALSE, 
                    useENSR=FALSE, genome=NULL, ...) { 

  # sort out assemblies
  g <- unique(genome(x))
  if (is.null(g)) g <- genome 
  if (!g %in% c("hg19","GRCh37","hg38","GRCh38")) {
    message("genome(x) is set to ",unique(genome(x)),", which is unsupported.")
    stop("Provide a `genome` argument, or set genome(x) manually, to proceed.")
  }
  if (!shrink) { 
    data("horvathAge", package="biscuiteer") 
    horvath <- horvathAge[[g]]
    intercept <- horvathAge$intercept 
  } else { 
    data("horvathAgeShrunken", package="biscuiteer") 
    horvath <- horvathAgeShrunken[[g]]
    intercept <- horvathAgeShrunken$intercept 
  } 
  seqlevelsStyle(horvath) <- seqlevelsStyle(x)

  # change padding if requested
  stopifnot(pad > 0)
  if (pad > 1) horvath <- trim(resize(horvath, pad*2, fix="center"))
  if (useENSR) {
    hasENSR <- which(!is.na(horvath$ENSR))
    message("Using ENSEMBL boundaries for ", length(hasENSR), " regions...") 
    start(horvath[hasENSR]) <- horvath[hasENSR]$ENSRstart
    end(horvath[hasENSR]) <- horvath[hasENSR]$ENSRend
  }

  # assess coverage (since this affects the precision of estimates)
  message("Assessing coverage across age-associated regions...") 
  covgWGBSage <- getCoverage(x, regions=horvath, what="perRegionTotal")
  rownames(covgWGBSage) <- horvath$name
  NAs <- which(covgWGBSage < minCovg, arr.ind=TRUE)

  if (nrow(NAs) > 0) {

    # warn the user if insufficient coverage is detected across a region(s)
    bySample <- split(horvath$name[NAs[,1]], sampleNames(x)[NAs[, 2]])
    message("Regions with coverage < ", minCovg, " in each sample:")
    pctMissing <- round((sapply(bySample, length)/length(horvath))*100)
    for (i in names(pctMissing)) message(i, ": ", pctMissing[i], "%")

    # warn the user if insufficient samples exist to impute across region(s) 
    subM <- DelayedMatrixStats::rowSums2(covgWGBSage>=minCovg,na.rm=T) < minSamp
    message(round(100*sum(subM)/length(subM)), "% of regions lack ",
            "sufficient coverage in enough samples to impute.")

    # is there a significant imbalance in positive/negative signes?
    tst <- fisher.test(table(horvath$sign, subM))
    message("Fisher's exact test (assessing sign bias due to missingness):")
    message("Odds ratio: ", round(tst$estimate,3), 
            " (p-value: ", round(tst$p.value,3), ") -- ",
            ifelse(tst$p.value < 0.1, "likely", "unlikely"), 
            " to significantly bias estimates.")

    # either way, probably a good idea to fix stuff 
    message("Raise `pad` (", pad, "), lower `minCovg` (", minCovg, "), ",
            "set `impute` and/or `useENSR` to TRUE.")
  } else { 
    message("All regions in all samples appear to be sufficiently covered.") 
  }

  # extract regional methylation estimates (and set < minCovg sites to NA) 
  methWGBSage <- getMeth(x, regions=horvath, type="raw", what="perRegion")
  rownames(methWGBSage) <- horvath$name
  methWGBSage[covgWGBSage < minCovg] <- NA
  methWGBSage <- as(methWGBSage, "matrix") # 353 x ncol(x) is not too huge

  # impute, if requested, any sites with insufficient coverage
  if (impute) methWGBSage <- impute.knn(methWGBSage, k=minSamp, ...)$data

  keep <- rowSums2(is.na(methWGBSage)) < 1
  tst <- fisher.test(table(horvath$sign, keep))
  message("Fisher's exact test (assessing sign bias due to dropped regions):")
  message("Odds ratio: ", round(tst$estimate,3), 
          " (p-value: ", round(tst$p.value,3), ") -- ",
          ifelse(tst$p.value < 0.1, "likely", "unlikely"), 
          " to significantly bias estimates.")
  methWGBSage <- methWGBSage[which(keep), ] 

  design <- rbind(Intercept=rep(1, ncol(x)), methWGBSage)
  coefs <- c(intercept, horvath[rownames(WGBSage)]$score)
  res <- list(age=(t(design) %*% coefs)[,1],
              coefs=coefs,
              meth=methWGBSage,
              call=sys.call())
  return(res)
    
}
