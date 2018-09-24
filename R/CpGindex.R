#' Measure hypermethylation-at-PRCs-in-CGIs and hypomethylation-at-WCGWs-in-PMDs
#'
#' At some point in some conference call somewhere, a collaborator suggested
#' that a simple index of Polycomb repressor complex (PRC) binding site hyper-
#' methylation and CpG-poor "partially methylated domain" (PMD) hypomethylation
#' would be a handy yardstick for both deterministic and stochastic changes 
#' associated with proliferation, aging, and cancer. This function provides 
#' such an index by compiling measures of aberrant hyper- and hypo-methylation
#' along with the ratio of hyper- to hypo-methylation. (The logic for this is 
#' that while the phenomena tend to occur together, there are many exceptions)
#' The resulting measures can provide a high-level summary of proliferation-,
#' aging-, and/or disease-associated changes in DNA methylation across samples.
#' 
#' The choice of defaults is fairly straightforward: in 2006, three independent 
#' groups reported recurrent hypermethylation in cancer at sites marked by both
#' H3K4me3 (activating) and H3K27me3 (repressive) histone marks in embryonic
#' stem cells; these became known as "bivalent" sites. The Roadmap Epigenome 
#' project performed ChIP-seq on hundreds of normal primary tissues and cell 
#' line results from the ENCODE project to generate a systematic catalog of 
#' "chromatin states" alongside dozens of whole-genome bisulfite sequencing 
#' experiments in the same tissues. We used both to generate a default atlas 
#' of bivalent (Polycomb-associated and transcriptionally-poised) sites from 
#' H9 human embryonic stem cells which retain low DNA methylation across normal
#' (non-placental) REMC tissues. In 2018, Zhou and Dinh (Nature Genetics) found
#' isolated [AT]CG[AT] sites, or "solo-WCGW" motifs, in common PMDs as the most
#' universal barometer of proliferation- and aging-associated methylation loss 
#' in mammalian cells, so we use their solo-WCGW sites in common PMDs as the 
#' default measure for hypomethylation. The resulting CpGindex is a vector of
#' length 3 for each sample: hypermethylation, hypomethylation, and their ratio.
#'
#' We suggest fitting a model for the composition of bulk samples (tumor/normal,
#' tissue1/tissue2, or whatever is most appropriate) prior to drawing any firm 
#' conclusions from the results of this function. For example, a mixture of 
#' two-thirds normal tissue and one-third tumor tissue may produce the same or 
#' lower degree of hyper/hypomethylation than high-tumor-content cell-free DNA 
#' samples from the blood plasma of the same patient. Intuition is simply not 
#' a reliable guide in such situations, which occur with some regularity. If 
#' orthogonal estimates of purity/composition are available (flow cytometry, 
#' ploidy, yield of filtered cfDNA), it is a Very Good Idea to include them. 
#' 
#' The default for this function is to use the HMM-defined CpG islands from 
#' Hao Wu's paper (Wu, Caffo, Jaffee, Irizarry & Feinberg, Biostatistics 2010) 
#' as generic "hypermethylation" targets inside of "bivalent" (H3K27me3+H3K4me3)
#' sites (identified in H9 embryonic stem cells & unmethylated across normals),
#' and the solo-WCGW sites within common partially methylated domains from 
#' Wanding Zhou and Huy Dinh's paper (Zhou, Dinh, et al, Nat Genetics 2018)
#' as genetic "hypomethylation" targets (as above, obvious caveats about tissue
#' specificity and user-supplied possibilities exist, but the defaults are sane
#' for many purposes, and can be exchanged for whatever targets a user wishes). 
#' 
#' The function returns all three components of the "CpG index", comprised of 
#' hyperCGI and hypoPMD (i.e. hyper, hypo, and their ratio). The PMD "score" is 
#' a base-coverage-weighted average of losses to solo-WCGW bases within PMDs;
#' the PRC score is similarly base-coverage-weighted but across HMM CGI CpGs,
#' within polycomb repressor complex sites (by default, the subset of state 23
#' segments in the 25-state, 12-mark ChromImpute model for H9 which have less 
#' than 10 percent CpG methylation across the CpG-island-overlapping segment in
#' all normal primary tissues and cells from the Reference Epigenome project). 
#' By providing different targets and/or regions, users can customize as needed.
#'
#' The return value is a CpGindex object, which is really just a DataFrame that 
#' knows about the regions at which it was summarized, and reminds the user of 
#' this when they implicitly call the `show` method on it.
#' 
#' @param x       a BSseq object
#' @param CGIs    a GRanges of CpG island regions, or NULL (default is HMM CGIs)
#' @param PRCs    a GRanges of Polycomb targets, or NULL (H9 state 23 low-meth)
#' @param WCGW    a GRanges of solo-WCGW sites, or NULL (default is PMD WCGWs)
#' @param PMDs    a GRanges of hypomethylating regions, or NULL (default PMDs) 
#'
#' @return        a CpGindex (DataFrame w/cols `hyper`, `hypo`, `ratio` + 2 GRs)
#' 
#' @import DelayedMatrixStats
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @import S4Vectors
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
  if (is.null(PRCs)) PRCs <- .fetch(x, "H9state23unmeth", suffix)
  hyperMeth <- .subsettedWithin(x, y=CGIs, z=PRCs) 

  # summarize hypomethylation (at WCGWs) by region
  message("Computing hypomethylation indices...") 
  if (is.null(PMDs)) PMDs <- .fetch(x, "PMDs", suffix)
  if (is.null(WCGW)) WCGW <- .fetch(x, "Zhou_solo_WCGW_inCommonPMDs", suffix)
  hypoMeth <- .subsettedWithin(x, y=WCGW, z=PMDs) 

  # summarize both via ratios
  message("Computing indices...") 
  res <- new("CpGindex", 
             hyperMethRegions=PRCs, 
             hypoMethRegions=PMDs,
             DataFrame(hyper=hyperMeth, 
                       hypo=hypoMeth, 
                       ratio=hyperMeth/hypoMeth))
  return(res)

}

# class definition
setClass("CpGindex", contains="DataFrame",
         slots=c(hyperMethRegions="GenomicRanges", 
                 hypoMethRegions="GenomicRanges"))

# default show method 
setMethod("show", "CpGindex", 
  function(object) {
    callNextMethod()
    if (length(slot(object, "hyperMethRegions")) > 0 | 
        length(slot(object, "hypoMethRegions")) > 0) {
      cat("  -------\n")
      cat("This object is just a DataFrame that",
          "has an idea of where it came from:\n")
    }
    if (length(slot(object, "hyperMethRegions")) > 0) { 
      cat("Hypermethylation was tallied across", 
          length(slot(object, "hyperMethRegions")), "regions (see", 
          paste0(as.character(match.call()[[2]]), "@hyperMethRegions)."), "\n")
    }
    if (length(slot(object, "hypoMethRegions")) > 0) {
      cat("Hypomethylation was tallied across", 
          length(slot(object, "hypoMethRegions")), "regions (see", 
          paste0(as.character(match.call()[[2]]), "@hypoMethRegions)."), "\n")
    }
  }) 

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
