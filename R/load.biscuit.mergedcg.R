#' a bsseq loader for Biscuit mergecg output (BED-like format, 3 cols/sample)
#' e.g. P01-028-T06.joint.cg.merged.hg19.bed.gz has 3 samples, and thus 12 cols
#'
#' @param filename    the file (compressed or not, doesn't matter) to load
#' @param sampleNames sample names for the bsseq object (if NULL, will create)
#' @param hdf5        make the object HDF5-backed? (FALSE) 
#' 
#' @return            a BSseq object from the bsseq package
#'
#' @import bsseq
#' @import data.table
#' @import HDF5Array
#'
#' @export
mergedCgToBsSeq <- function(filename, sampleNames=NULL, hdf5=FALSE) {

  input <- filename
  if (base::grepl(".gz$", filename)) input <- paste("zcat", input)
  if (base::grepl(".bz2$", filename)) input <- paste("bzcat", input)

  if (!is.null(sampleNames)) {
    nsamples <- length(sampleNames)
    ncolumns <- 3 + (3*nsamples)
    dropcols <- seq(6, ncolumns, 3)
    merged.dt <- fread(input, sep="\t", sep2=",", na.string=".", drop=dropcols)
  } else {
    merged.dt <- fread(input, sep="\t", sep2=",", na.string=".") 
    nsamples <- (ncol(merged.dt) - 3) / 3
    ncolumns <- 3 + (3*nsamples)
    dropcols <- seq(6, ncolumns, 3)
    kept <- setdiff(seq_len(ncolumns), dropcols)
    merged.dt <- merged.dt[, kept, with=FALSE]
    sampleNames <- paste0("sample", seq_len(nsamples))
  }

  cols <- c("chr","start","end")
  sampcols <- apply(expand.grid(c("beta","covg"), seq_len(nsamples)), 
                    1, paste, collapse="")
  colnames(merged.dt) <- c(cols, sampcols)
  merged.dt[, "start"] <- merged.dt[, "start"] + 1 # quirk
  betacols <- paste0("beta", seq_len(nsamples))
  covgcols <- paste0("covg", seq_len(nsamples))
  message("Loaded data from ", filename, ". Creating bsseq object...")

  if (hdf5) { 
    hdf5_M <- writeHDF5Array(fixNAs(round(merged.dt[, betacols, with=FALSE]* 
                                          merged.dt[, covgcols, with=FALSE])))
    hdf5_Cov <- writeHDF5Array(fixNAs(merged.dt[, covgcols, with=FALSE]))
    BSseq(gr=makeGRangesFromDataFrame(merged.dt[, c("chr","start","end")]),
          M=hdf5_M, Cov=hdf5_Cov, sampleNames=sampleNames, rmZeroCov=TRUE)
  } else { 
    BSseq(gr=makeGRangesFromDataFrame(merged.dt[, c("chr","start","end")]),
          M=fixNAs(round(merged.dt[, betacols, with=FALSE]* 
                         merged.dt[, covgcols, with=FALSE])),
          Cov=fixNAs(merged.dt[, covgcols, with=FALSE]), 
          sampleNames=sampleNames, 
          rmZeroCov=TRUE)
  } 
}
