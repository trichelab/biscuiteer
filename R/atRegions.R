#' Summarize a bsseq dataset over defined regions
#'
#' Calls summarizeBsSeqOver to summarize a bsseq object over provided DNA
#' regions. Useful for exploring genomic data using cBioPortal.
#'
#' @param bsseq     A bsseq object
#' @param regions   A GRanges or GRangesList of regions
#' @param mappings  A mapping table with rownames(mappings) == colnames(bsseq)
#'                    (DEFAULT: NULL)
#' @param nm        Column of the mapping table to map to
#'                    (DEFAULT: "POETICname")
#' @param ...       Other arguments to pass to summarizeBsSeqOver
#'
#' @return          Summarized information about the bsseq object for the given
#'                    DNA regions
#'
#' @examples
#'
#' @export
#'
atRegions <- function(bsseq,
                      regions,
                      mappings = NULL,
                      nm = "POETICname",
                      ...) { 
  regional <- summarizeBsSeqOver(bsseq, regions, ...)
  regions <- regions[which(as.character(regions) %in% rownames(regional))]
  regional <- regional[as.character(regions), ]
  rownames(regional) <- names(regions) 
  if (!is.null(mappings)) colnames(regional) <- mappings[colnames(regional), nm]
  return(regional)
}
