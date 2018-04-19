#' a bsseq loader for Biscuit cg/ch output (BED-like format, 2 cols/sample)
#' e.g. P01-028-T06.joint.ch.hg19.bed.gz has 3 samples, and thus 9 columns
#'
#' @param filename    the file (compressed or not, doesn't matter) to load
#' @param sampleNames sample names (if NULL, create; if data.frame, make pData)
#' @param hdf5        make the object HDF5-backed? (FALSE; use in-core storage) 
#' @param sparse      make the object Matrix-backed? (TRUE)
#' 
#' @return            a BSseq object from the bsseq package
#'
#' @import bsseq
#' @import data.table
#' @import GenomicRanges
#'
#' @concept bsseq
#' @concept BSseq 
#'
#' @seealso load.biscuit
#' @seealso checkBiscuitBED
#'
#' @export
load.biscuit.unmerged <- function(filename,
                                  sampleNames=NULL,
                                  hdf5=FALSE,
                                  sparse=FALSE) {

  params <- checkBiscuitBED(filename, sampleNames, merged=FALSE)
  chh <- ifelse(base::grepl("c(p?)g", filename, ignore=TRUE), "CpG", "CpH")
  filename <- params$tbx$path
  message("Reading unmerged ", chh, " input from ", filename, "...")

  if (params$passes > 1) {
    # for files that tend to fill up /tmp or /shm, scan as tabix, yielding
    to.dt <- function(elt) data.table(read.table(textConnection(elt), sep="\t"))
    unmerged.dt <- data.table(params$preamble)[0,]
    loci <- 0
    while(length(res <- scanTabix(params$tbx)[[1]])) {
      loci <- loci + length(res)
      unmerged.dt <- rbind(unmerged.dt, Map(to.dt, res))
      message(loci, " ", chh, " loci processed")
    }
  } else { 
    unmerged.dt <- fread(paste("zcat", params$tbx$path), sep="\t", sep2=",", 
                         na.string=".", skip=ifelse(hasHeader, 1, 0))
  }
  colnames(unmerged.dt) <- params$colNames
  unmerged.dt[, "start"] <- unmerged.dt[, "start"] + 1 # quirk
  message("Loaded data from ", filename, ". Creating bsseq object...")

  if (hdf5) { 
    with(params, 
         makeBSseq_hdf5(unmerged.dt, 
                        betacols, 
                        covgcols, 
                        pData=pData, 
                        sparse=sparse))
  } else { 
    with(params, 
         makeBSseq(unmerged.dt, 
                   betacols, 
                   covgcols, 
                   pData=pData,
                   sparse=sparse))
  } 

}
