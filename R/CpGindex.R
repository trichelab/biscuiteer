#' simple index of hypermethylation-at-CpG-islands vs hypomethylation-in-WCGWs
#'
#' The default for this function is to use the HMM-defined CpG islands from 
#' Hao Wu's paper (Wu, Caffo, Jaffee, Irizarry & Feinberg, Biostatistics 2010) 
#' as generic "hypermethylation" targets (obviously tissue-specific subsets 
#' might be more informative on a number of levels, e.g. as bivalent sites),
#' and the solo-WCGW sites within common partially methylated domains from 
#' Wanding Zhou and Huy Dinh's paper (Zhou, Dinh, et al, Nat Genetics 2018)
#' as genetic "hypomethylation" targets (as above, obvious caveats about tissue
#' specificity and user-supplied possibilities exist). 
#' 
#' The return value of the function is an idiotically simple measure, comprised
#' of hyperCGI/hypoPMD (both the components and the ratio). The PMD "score" is 
#' a base-coverage-weighted average of losses to solo-WCGW bases within PMDs;
#' the CGI score is similarly base-coverage-weighted but across HMM CGI CpGs. 
#' 
#' By providing different targets a user can customize as needed.
#' 
#' @param x       a BSseq object
#' @param CGIs    a GRanges of hypermethylation targets, or NULL (default CGIs)
#' @param WCGW    a GRanges of hypomethylation targets, or NULL (default WCGW)
#'
#' @return        a list: `hyper`, `hypo`, and `ratio` (data.frame) by sample
#' 
#' @import GenomicRanges
#' @import GenomeInfoDb
#' 
#' @export 
CpGindex <- function(x, CGIs=NULL, WCGW=NULL, PMDs=NULL) {

  if (is.null(unique(genome(x)))) { 
    stop("You must assign a genome to your BSseq object")
  } else { 
    genome <- unique(genome(x))
    if (genome %in% c("hg19","GRCh37")) suffix <- "hg19"
    else if (genome %in% c("hg38","GRCh38")) suffix <- "hg38"
    else stop("Only human genomes (hg19/GRCh37, hg38/GRCh38) are supported ATM")
  } 

  if (is.null(CGIs)) {
    dat <- paste0("HMM_CpG_islands.", suffix)
    message("Loading ", dat, "...") 
    data(list=dat, package="biscuiteer") 
    CGIs <- get(dat)
    seqlevelsStyle(CGIs) <- seqlevelsStyle(x)
  }
  message("Computing hypermethylation indices...") 
  hyper <- getMeth(x, regions=CGIs, type="raw", what="perRegion")

  if (is.null(WCGW)) { 
    dat <- paste0("Zhou_solo_WCGW_inCommonPMDs.", suffix)
    message("Loading ", dat, "...") 
    data(list=dat, package="biscuiteer") 
    WCGW <- get(dat) 
    seqlevelsStyle(WCGW) <- seqlevelsStyle(x)
  }
  message("Computing hypomethylation indices...") 
  hypo <- getMeth(x, regions=WCGW, type="raw", what="perRegion")

  message("Computing ratio...") 
  ratio <- DataFrame(hyper=DelayedArray::colMeans(hyper, na.rm=TRUE),
                     hypo=DelayedArray::colMeans(hypo, na.rm=TRUE))
  ratio$ratio <- ratio$hyper / ratio$hypo

  res <- list(hyper=hyper, hypo=hypo, ratio=ratio)
  return(res) 

} 
