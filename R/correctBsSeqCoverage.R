#' Tweaked version of TitanCNA's preprocessing
#'
#' Limits coverage to chromosomes 1 through 22 and X (for better or for worse)
#'
#' @param tumor   Binned tumor coverage - a GRanges from binCoverage
#' @param normal  Binned normal coverage - a GRanges from binCoverage
#' @param bsgc    Binned c2t GC content (DEFAULT: bsgc.hg19, 1 kb bins)
#' @param maps    Binned c2t mappability (DEFAULT: maps.hg19, 1 kb bins)
#'
#' @return        Corrected tumor and normal read counts
#'
#' @importFrom utils data
#' @importFrom GenomeInfoDb keepSeqlevels
#' @import GenomicRanges
#' @import HMMcopy
#'
#' @seealso binCoverage
#' @seealso TitanCNA::correctReadDepth
#'
#' @examples
#'
#'   tum_bed <- system.file("extdata", "MCF7_Cunha_chr11p15_shuffled.bed.gz",
#'                          package="biscuiteer")
#'   nor_bed <- system.file("extdata", "MCF7_Cunha_chr11p15.bed.gz",
#'                          package="biscuiteer")
#'   tum_vcf <- system.file("extdata", "MCF7_Cunha_shuffled_header_only.vcf.gz",
#'                          package="biscuiteer")
#'   nor_vcf <- system.file("extdata", "MCF7_Cunha_header_only.vcf.gz",
#'                          package="biscuiteer")
#'
#'   tumor <- read.biscuit(BEDfile=tum_bed, VCFfile=tum_vcf, merged=FALSE)
#'   normal <- read.biscuit(BEDfile=nor_bed, VCFfile=nor_vcf, merged=FALSE)
#'
#'   corr <- correctBsSeqCoverage(tumor=tumor, normal=normal)
#'
#' @export
#'
correctBsSeqCoverage <- function(tumor,
                                 normal,
                                 bsgc = NULL,
                                 maps = NULL) {

  if (is.null(bsgc)) { # FIXME: support hg38
    data(bsgc.hg19) 
    bsgc <- bsgc.hg19
  }
  if (is.null(maps)) { # FIXME: support hg38
    data(maps.hg19) 
    maps <- maps.hg19
  }
  if (is.null(attr(tumor, "binned"))) tumor <- binCoverage(tumor, bsgc)
  if (is.null(attr(normal, "binned"))) normal <- binCoverage(normal, bsgc)

  normal <- subsetByOverlaps(normal, tumor)
  tumor <- subsetByOverlaps(tumor, normal)
  stopifnot(length(tumor) == length(normal))

  bsgc <- subsetByOverlaps(bsgc, tumor)
  stopifnot(length(tumor) == length(bsgc))
  maps <- subsetByOverlaps(maps, tumor)
  stopifnot(length(maps) == length(bsgc))
  tumor$gc <- normal$gc <- bsgc$score
  tumor$map <- normal$map <- maps$score
  names(mcols(tumor)) <- c("reads","gc","map")
  names(mcols(normal)) <- c("reads","gc","map")
  
  message("Correcting tumor for bsgc bias and mappability...")
  tumor_copy <- correctDepth(tumor)
    
  message("Normalizing tumor against normal sample...")
  tumor_copy$copy <- tumor_copy$copy - correctDepth(normal)$copy
  tumor_copy$score <- tumor_copy$copy
  seqlevels(tumor_copy) <- seqlevels(bsgc) # Fill in missing seqlevels
  tumor_copy <- keepSeqlevels(tumor_copy, 
                              paste0("chr", c(1:22, "X")), 
                              pruning.mode="coarse")
  return(tumor_copy) 
}
