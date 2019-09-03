#' Matched tumor-normal coverage processing and segmentation
#'
#' @param stub0  bigWig base file name for normal (i.e. stub0.covg.hg19.bw)
#' @param stub1  bigWig base file name for tumor (i.e. stub1.covg.hg19.bw)
#'                 (DEFAULT: NULL)
#' @param bsgc   Binned c2t GC content (DEFAULT: bsgc.hg19, 1 kb bins)
#' @param maps   Binned c2t mappability (DEFAULT: maps.hg19, 1 kb bins)
#' @param ...    Other parameters to pass to WGBSseg
#'
#' @return       A GRanges object
#'
#' @import fastseg
#' @import Biostrings
#'
#' @seealso WGBSseg
#'
#' @examples
#'
#'   stubs <- c("P01-010", "P01-012", "P01-015", "P01-016", "P01-017")
#'   stubs <- c(stubs, "P01-018", "P01-020", "P01-021")
#'   names(stubs) <- paste0(stubs, "-T06-T-P01")
#'   names(stubs)[1] <- "P01-010-T07-T-P01"
#'   segs <- list()
#'   for(i in names(stubs)) segs[[i]] <- TN_WGBS_CNA(stubs[[i]])$segs # robust
#'   segs <- GRangesList(segs) # to feed to grToSeg and get full .seg output
#'   for(i in names(segs)) export(segs[[i]], paste0(i, ".CN.hg19.bw"))
#'
#' @export
#'
TN_WGBS_CNA <- function(stub0,
                        stub1 = NULL,
                        bsgc = NULL,
                        maps = NULL,
                        ...) {
  if (is.null(stub1)) {
    stub <- stub0
    stub0 <- paste0(stub, "-normal")
    stub1 <- paste0(stub, "-tumor")
  } else { 
    stubs <- c(stub0, stub1)
    stub <- strsplit(consensusString(stubs, thresh=1), '-?', fixed=TRUE)[[1]][1]
  }
  file0 <- paste0(stub0, ".covg.hg19.bw")
  file1 <- paste0(stub1, ".covg.hg19.bw")
  message("\n", "Processing coverage files for ", stub, "...", "\n")
  files <- c(normal=file0, tumor=file1)
  if (all(c("tumor","normal") %in% names(files))) {
    message("Reading ", paste(files, collapse=" and "), " ...")
    covgs <- GRangesList(sapply(files[c("tumor","normal")], import))
    if (seqlengths(covgs)["chrX"] != 155270560) {
      stop("These don't look like hg19 tracks. Aborting!")
    } else {
      cn <- correctBsSeqCoverage(covgs$tumor, covgs$normal, bsgc, maps)
      message("Segmenting...")
      segs <- WGBSseg(cn, ...)
      message("Done.") 
      res <- list(cn=cn, segs=segs)
    }
  } else {
    message("Could not find matched tumor-normal pair for ", stub)
    message("Found the following files matching ", stub, ":")
    for (i in names(files)) {
      message(i, ": ", files[i])
    }
    stop("Check your files and try again, or run WGBSseg manually.")
  }
  return(res)
}
