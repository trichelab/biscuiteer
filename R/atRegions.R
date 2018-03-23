#' summarize a bsseq dataset over defined regions (for e.g. cBioPortal) 
#' 
#' @param bsseq     a bsseq object
#' @param regions   a GRanges or GRangesList of regions 
#' @param mappings  a mapping table, with rownames(mappings) == colnames(bsseq)
#' @param nm        column of the mapping table to map to (default "POETICname")
#'
#' @import bsseq
#'
#' @export
atRegions <- function(bsseq, regions, mappings=NULL, nm="POETICname") { 
  regional <- summarizeBsSeqOver(bsseq, regions)[as.character(regions), ]
  regions <- regions[which(as.character(regions) %in% rownames(regional))]
  rownames(regional) <- names(regions) 
  if (!is.null(mappings)) colnames(regional) <- mappings[colnames(regional), nm]
  return(regional)
}
