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
#'   tcga_bed <- system.file("extdata", "TCGA_BLCA_A13J_chr11p15_merged.bed.gz",
#'                           package = "biscuiteer")
#'   tcga_vcf <- system.file("extdata", "TCGA_BLCA_A13J_header_only.vcf.gz",
#'                           package = "biscuiteer")
#'   bisc <- read.biscuit(BEDfile = tcga_bed, VCFfile = tcga_vcf,
#'                        merged = TRUE, genome = "hg38", verbose = TRUE)
#'
#'   reg <- GRanges(seqnames = rep("chr11",5),
#'                  strand = rep("*",5),
#'                  ranges = IRanges(start = c(0,2.8e6,1.17e7,1.38e7,1.69e7),
#'                                   end= c(2.8e6,1.17e7,1.38e7,1.69e7,2.2e7))
#'                  )
#'   regions <- atRegions(bsseq = bisc, regions = reg)
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
