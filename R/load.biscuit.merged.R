#' a bsseq loader for Biscuit mergecg output (BED-like format, 3 cols/sample)
#' e.g. P01-028-T06.joint.cg.merged.hg19.bed.gz has 3 samples, and thus 12 cols
#'
#' @param params      parameters from having run checkBiscuitBED 
#' 
#' @return            a BSseq object from the bsseq package
#'
#' @import bsseq
#' @import readr
#'
#' @seealso load.biscuit.unmerged
#'
#' @export
load.biscuit.merged <- function(params) {

  ncolumns <- 3 + (3 * params$nSamples)
  dropcols <- seq(6, ncolumns, 3)
  filename <- params$tbx$path
  tblCols <- params$colNames[setdiff(seq_along(params$colNames), dropcols)]
  message("Reading merged CpG input from ", filename, "...")
  tbl <- read_tsv(filename, na.string=".", comment="#", 
                     drop=dropcols, skip=ifelse(params$hasHeader, 1, 0))
  colnames(merged.dt) <- dtCols
  merged.dt[, "start"] <- merged.dt[, "start"] + 1 # quirk
  message("Loaded data from ", filename, ". Creating bsseq object...")

  if (params$hdf5) { 
    with(params, makeBSseq_hdf5(merged.dt, betacols, covgcols, pData, sparse))
  } else { 
    with(params, makeBSseq(merged.dt, betacols, covgcols, pData, sparse))
  } 
}
