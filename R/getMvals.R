#' helper function for compartment inference (shrink-by-smoothing)
#' 
#' since we want something resembling M-values for compartment inference,
#' this function grabs suitably deep (>= 3 reads in >=1 sample) measurements
#' and turns them into lightly moderated M-values for compartment calling
#' by Dirichlet smoothing (adding addReads reads to each of M and Covg)
#' 
#' @param bsseq       a BSseq object with methylated and total reads 
#' @param chrom       what chromosome to grab ("chr22")
#' @param minCov      minimum read coverage to use 
#' @param resolution  resolution for binning (1e5)
#' @param k           pseudoreads for smoothing
#' 
#' @return          a matrix of Mvalues
#'
#' @import bsseq
#'
#' @export
getMvals <- function(bsseq, chrom="chr22", minCov=3, resolution=1e5, k=1) { 

  usable <- seqnames(bsseq) == chrom 
  usable[which(usable)] <- rowSums(getCoverage(bsseq[usable,]) >= minCov) > 0
  SummarizedExperiment(
    assays=list(Mval=flogit(fixNAs(getSmoothedBeta(bsseq[usable,], k=k)))), 
    rowRanges=granges(bsseq)[usable]
  )

}

# helper fn
getSmoothedBeta <- function(x, k=1) (assays(x)$M + k)/(assays(x)$Cov + (k * 2))
