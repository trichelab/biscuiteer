#' Helper function for compartment inference
#'
#' Want an object with nominally Gaussian error for compartment inference, so
#' this function uses 'suitable' (defaults to to 3 or more reads in 2 or more 
#' samples) measurements. Using Dirichlet smoothing (adding 'k' reads to M
#' and U), these measurements are then turned into lightly moderated,
#' logit-transformed methylated-fraction estimates for compartment calling.
#'
#' @param x        A bsseq object with methylated and total reads
#' @param minCov   Minimum read coverage for landmarking samples (DEFAULT: 3)
#' @param minSamp  Minimum landmark samples with at least minCov (DEFAULT: 2)
#' @param k        Pseudoreads for smoothing (DEFAULT: 0.1)
#' @param r        Regions to collapse over - if NULL, do it by CpG
#'                   (DEFAULT: NULL)
#'
#' @return         Smoothed logit(M / Cov) matrix with coordinates as row names
#'
#' @import gtools
#' @import bsseq
#'
#' @aliases getMvals
#'
#' @examples
#'
#'   tcga_bed <- system.file("extdata", "TCGA_BLCA_A13J_chr11p15_merged.bed.gz",
#'                           package = "biscuiteer")
#'   tcga_vcf <- system.file("extdata", "TCGA_BLCA_A13J_header_only.vcf.gz",
#'                           package = "biscuiteer")
#'   bisc <- read.biscuit(BEDfile = tcga_bed, VCFfile = tcga_vcf,
#'                        merged = TRUE, genome = "hg38")
#'
#'   reg <- GRanges(seqnames = rep("chr11",5),
#'                  strand = rep("*",5),
#'                  ranges = IRanges(start = c(0,2.8e6,1.17e7,1.38e7,1.69e7),
#'                                   end= c(2.8e6,1.17e7,1.38e7,1.69e7,2.2e7))
#'                  )
#'
#'   frac <- getLogitFracMeth(bisc, minSamp = 1, r = reg)
#'
#' @export
#'
getLogitFracMeth <- function(x,
                             minCov = 3,
                             minSamp = 2,
                             k = 0.1,
                             r = NULL) {

  # do any loci/regions have enough read coverage in enough samples? 
  if (!is.null(r) && is(r, "GenomicRanges")) {
    covgs <- getCoverage(x, sort(r), type="Cov", what="perRegionTotal")
  } else { 
    covgs <- getCoverage(x, type="Cov", what="perBase")
  } 

  usable <- DelayedMatrixStats::rowSums2(covgs >= minCov) >= minSamp
  if (!any(usable)) stop("No usable loci/regions ( >= minCov in >= minSamp )!")
    
  # construct a subset of the overall BSseq object with smoothed mvalues 
  if (!is.null(r) && is(r, "GenomicRanges")) {
    getSmoothedLogitFrac(x, k=k, minCov=minCov, r=subset(sort(r), usable))
  } else { 
    getSmoothedLogitFrac(subset(x, usable), k=k, minCov=minCov)
  } 

}

# Helper function to find smoothed logit fraction
# x: bsseq object
# k: pseudoreads for smoothing (DEFAULT: 0.1)
# minCov: minimum coverage (DEFAULT: 3)
# maxFrac: maximum fraction of NAs allowed (DEFAULT 0.5)
# r: regions to collapse over
getSmoothedLogitFrac <- function(x,
                                 k = 0.1,
                                 minCov = 3,
                                 maxFrac = 0.5,
                                 r = NULL) {

  if (!is.null(r) && is(r, "GenomicRanges")) {
    M <- getCoverage(x, sort(r), type="M", what="perRegionTotal")
    U <- getCoverage(x, sort(r), type="Cov", what="perRegionTotal") - M 
    rnames <- as.character(sort(r))
  } else { 
    M <- getCoverage(x, type="M", what="perBase")
    U <- getCoverage(x, type="Cov", what="perBase") - M 
    rnames <- as.character(granges(x))
  } 

  res <- logit((M + k) / ((M + k) + (U + k))) 
  rownames(res) <- rnames 

  makeNA <- ((M + U) < minCov)
  maxPct <- paste0(100 * maxFrac, "%")
  tooManyNAs <- (DelayedMatrixStats::colSums2(makeNA)/nrow(x)) > maxFrac
  if (any(tooManyNAs)) {
    message(paste(colnames(x)[tooManyNAs],collapse=", ")," are >",maxPct," NA!")
  }
  res[ makeNA ] <- NA
  return(res)

}

#' @describeIn getLogitFracMeth Alias for getLogitFracMeth
#'
#' @export
#'
getMvals <- getLogitFracMeth
