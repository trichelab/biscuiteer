#' A BED checker for Biscuit cg/ch output (BED-like format, 2 or 3 cols/sample).
#' By default, files with over 50 million loci will be processed iteratively,
#' since data.table tends to run into problems with .gzipped joint CpH files. 
#'
#' @param filename    the file (compressed or not, doesn't matter) to load
#' @param sampleNames sample names (if NULL, create; if data.frame, make pData)
#' @param merged      is the file a merged CpG file? (if NULL, will guess) 
#' @param sparse      make the object Matrix-backed? (NULL; do so if beneficial)
#' @param tabix       for files > `tabix` lines long, use TabixFile (5e7)
#' 
#' @return            parameters for makeBSseq or makeBSseq_hdf5
#'
#' @import LaF
#'
#' @seealso load.biscuit
#'
#' @export
checkBiscuitBED <- function(filename,
                            sampleNames=NULL,
                            merged=NULL, 
                            sparse=NULL,
                            tabix=5e7){

  input <- filename
  if (base::grepl(".gz$", filename)) input <- paste("zcat", input)
  if (base::grepl(".bz2$", filename)) input <- paste("bzcat", input)
  nlines <- determine_nlines(filename)
  use_tabix <- (nlines > tabix)

  # read the first few samples and see if we have any major problems
  preamble <- read.table(filename, sep="\t", na.strings=".", nrows=10)
  if (is.null(merged)) merged <- base::grepl("merged", ignore=TRUE, filename)
  colsPerSample <- ifelse(merged, 3, 2)
  nSamples <- (ncol(preamble) - 3) / colsPerSample
  if (!is.null(sampleNames)) {
    if (is(sampleNames, "DataFrame") | is(sampleNames, "data.frame")) {
      stopifnot(ncol(sampleNames) == nSamples)
      pData <- DataFrame(sampleNames)
    } else {
      stopifnot(length(sampleNames) == nSamples)
      pData <- DataFrame(sampleName=sampleNames)
    }
  } else {
    sampleNames <- paste0("sample", seq_len(nSamples))
    pData <- DataFrame(sampleName=sampleNames)
  }

  cols <- c("chr","start","end")
  sampcols <- apply(expand.grid(c("beta","covg"), seq_len(nSamples)), 
                    1, paste, collapse="")
  colNames <- base::gsub(" ", "", c(cols, sampcols)) # quirk
  betacols <- paste0("beta", seq_len(nSamples))
  covgcols <- paste0("covg", seq_len(nSamples))
  message(filename, " looks valid for import.")
  colnames(preamble) <- colNames

  params <- list(input=input,
                 preamble=preamble,
                 nSamples=nSamples,
                 colNames=colNames,
                 betacols=betacols,
                 covgcols=covgcols,
                 use_tabix=use_tabix,
                 nlines=nlines,
                 pData=pData,
                 yield=1e5)
  return(params) 

}
