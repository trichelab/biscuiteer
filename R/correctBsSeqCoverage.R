#' tweaked version of TitanCNA's preprocessing
#' 
#' limits coverage to chr{1:22,X}, for better or for worse
#'
#' @param tumor     binned tumor coverage (a GRanges from binCoverage)
#' @param normal    binned normal coverage (a GRanges from binCoverage)
#' @param bsgc      binned c2t GC content (default is bsgc.hg19, 1kb bins)
#' @param maps      binned c2t mappability (default is maps.hg19, 1kb bins)
#'
#' @return          corrected tumor and normal read counts
#'
#' @import GenomicRanges
#' @import HMMcopy
#' 
#' @export
correctBsSeqCoverage <- function(tumor, normal, bsgc=NULL, maps=NULL) {

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
  if (length(tumor) != length(bsgc)) browser()
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
  tumor_copy <- keepSeqlevels(tumor_copy, 
                              paste0("chr", c(1:22, "X")), 
                              pruning.mode="coarse")
  return(tumor_copy) 
}
