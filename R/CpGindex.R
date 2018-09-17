#' simple index of hypermethylation-at-CpG-islands vs hypomethylation-in-WCGWs
#'
#' The default for this function is to use the HMM-defined CpG islands from 
#' Hao Wu's paper (Wu, Caffo, Jaffee, Irizarry & Feinberg, Biostatistics 2010) 
#' as generic "hypermethylation" targets (obviously tissue-specific subsets 
#' might be more informative on a number of levels, e.g. PRC bivalent sites),
#' and the solo-WCGW sites within common partially methylated domains from 
#' Wanding Zhou and Huy Dinh's paper (Zhou, Dinh, et al, Nat Genetics 2018)
#' as genetic "hypomethylation" targets (as above, obvious caveats about tissue
#' specificity and user-supplied possibilities exist). 
#' 
#' The return value of the function is an idiotically simple measure, comprised
#' of hyperCGI/hypoPMD (both components, and the ratio). The PMD "score" is 
#' a base-coverage-weighted average of losses to solo-WCGW bases within PMDs;
#' the CGI score is similarly base-coverage-weighted but across HMM CGI CpGs,
#' within polycomb repressor complex sites if provided (by default, hESC PRCs; 
#' specifically, state 23 in the REMC imputed 25-state, 12-mark model for H9).
#' 
#' By providing different targets and/or regions, users can customize as needed.
#' 
#' @param x       a BSseq object
#' @param CGIs    a GRanges of CpG island regions, or NULL (default is HMM CGIs)
#' @param PRCs    a GRanges of Polycomb targets, or NULL (default H9 bivalent)
#' @param WCGW    a GRanges of solo-WCGW sites, or NULL (default is PMD WCGWs)
#' @param PMDs    a GRanges of hypomethylating regions, or NULL (default PMDs) 
#'
#' @return        a data.frame with columns `hyper`, `hypo`, and `ratio`
#' 
#' @import DelayedMatrixStats
#' @import GenomicRanges
#' @import GenomeInfoDb
#' 
#' @export 
CpGindex <- function(x, CGIs=NULL, PRCs=NULL, WCGW=NULL, PMDs=NULL) {

  # necessary evil 
  if (is.null(unique(genome(x)))) { 
    stop("You must assign a genome to your BSseq object before proceeding.")
  } else { 
    genome <- unique(genome(x))
    if (genome %in% c("hg19","GRCh37")) suffix <- "hg19"
    else if (genome %in% c("hg38","GRCh38")) suffix <- "hg38"
    else stop("Only human genomes (hg19/GRCh37, hg38/GRCh38) are supported ATM")
  } 

  # summarize hypermethylation by region 
  message("Computing hypermethylation indices...") 
  if (is.null(CGIs)) CGIs <- .fetch(x, "HMM_CpG_islands", suffix) 
  if (is.null(PRCs)) PRCs <- .fetch(x, "H9bivalent", suffix)
  hyperMeth <- .subsettedWithin(x, y=CGIs, z=PRCs) 

  # summarize hypomethylation (at WCGWs) by region
  message("Computing hypomethylation indices...") 
  if (is.null(PMDs)) PMDs <- .fetch(x, "PMDs", suffix)
  if (is.null(WCGW)) WCGW <- .fetch(x, "Zhou_solo_WCGW_inCommonPMDs", suffix)
  hypoMeth <- .subsettedWithin(x, y=WCGW, z=PMDs) 

  # summarize both via ratios
  message("Computing ratios...") 
  res <- data.frame(hyper=hyperMeth, hypo=hypoMeth, ratio=hyperMeth/hypoMeth)
  attr(res, "hyperMethRegions") <- PRCs
  attr(res, "hypoMethRegions") <- PMDs
  return(res)

}

# helper fn
.fetch <- function(x, prefix, suffix) {
  dat <- paste(prefix, suffix, sep=".")
  message("Loading ", dat, "...") 
  data(list=dat, package="biscuiteer") 
  xx <- get(dat)
  seqlevelsStyle(xx) <- seqlevelsStyle(x)
  return(xx)
} 

# helper fn
.subsettedWithin <- function(x, y, z) { 
  y <- sort(subsetByOverlaps(y, x))
  z <- sort(subsetByOverlaps(z, y))
  res <- getMeth(subsetByOverlaps(x,y), regions=z, type="raw", what="perRegion")
  if (any(is.na(res)) | any(is.nan(res))) {
    return(.delayedNoNaN(res, z)) # slower but exact 
  } else { 
    return(res %*% width(z)/sum(width(z))) # much faster
  }
}

# helper fn
.delayedNoNaN <- function(x, z) {
  res <- c() 
  for (i in colnames(x)) {
    use <- !is.na(x[,i]) & !is.nan(x[,i])
    zz <- z[which(as(use, "logical"))]
    res[i] <- as.matrix(x[use, i]) %*% (width(zz)/sum(width(zz)))
  }
  return(res)
}
