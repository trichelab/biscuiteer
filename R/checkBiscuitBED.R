#' A BED checker for Biscuit cg/ch output (BED-like format, 2 or 3 cols/sample).
#' By default, files with over 50 million loci will be processed iteratively,
#' since data.table tends to run into problems with .gzipped joint CpH files. 
#' This function absolutely assumes that BED files are tabixed. No exceptions!
#' 
#' @param filename    the file (compressed and tabixed, with header) to load
#' @param sampleNames sample names (if NULL, create; if data.frame, make pData)
#' @param yieldSize   for files > `yieldSize` lines long, chunk the file (5e7)
#' 
#' @return            parameters for makeBSseq or makeBSseq_hdf5
#'
#' @import            Rsamtools 
#'
#' @seealso           load.biscuit
#'
#' @export
checkBiscuitBED <- function(filename, 
                            sampleNames=NULL, 
                            yieldSize=5e7, 
                            merged=NULL){

  if (!base::grepl(".gz$", filename)) stop("Only tabix'ed files are supported.")
  message("Checking ", filename, " for import...")
  tbx <- TabixFile(filename, yieldSize=yieldSize)

  # look for Tabix header, then look for problems
  hasHeader <- (length(headerTabix(tbx)$header) > 0)
  if (hasHeader) { 

    # {{{ has a header; easy and verifiable for ASM
    message(filename," has a header line and is ", appendLF=FALSE)
    preamble <- read.table(tbx$path, header=TRUE, sep="\t", na.strings=".", 
                           comment.char="`", nrows=3)
    colNames <- base::sub("#", "", colnames(preamble))
    if (is.null(merged)) merged <- any(grepl("context", colnames(preamble)))
    message(ifelse(merged, "merged", "unmerged"), " data.") 

    betacols <- grep("beta", colnames(preamble), value=TRUE) 
    covgcols <- grep("covg", colnames(preamble), value=TRUE) 
    sNames <- condenseSampleNames(tbx, stride=ifelse(merged, 3, 2))
    if (is.null(sampleNames)) sampleNames <- sNames 
    nSamples <- length(sNames)
    # }}}

  } else {

    # {{{ no header; guesswork
    message(paste0(filename," has no header (!!!) and is "), appendLF=FALSE)
    preamble <- read.table(tbx$path, sep="\t", na.strings=".", nrows=3)
    if (is.null(merged)) merged <- grepl("merged", filename, ignore=TRUE)
    message(ifelse(merged, "merged", "unmerged"), " data.") 
    cols <- c("chr","start","end")
    colsPerSample <- ifelse(merged, 3, 2)
    nSamples <- (ncol(preamble) - 3) / colsPerSample
    sampleses <- paste0("sample", seq_len(nSamples))
    if (merged) colSuffixes <- c(".beta",".covg",".context")
    else colSuffixes <- c(".beta",".covg")
    sampcols <- paste0(rep(sampleses, each=ifelse(merged, 2, 3)), 
                       rep(colSuffixes, nSamples))
    colNames <- base::gsub(" ", "", c(cols, sampcols)) # quirk
    # }}}

  }

  if (!is.null(sampleNames)) { 
    
    # {{{ if sampleNames is provided OR deduced
    if (is(sampleNames, "DataFrame") | is(sampleNames, "data.frame")) {
      stopifnot(ncol(sampleNames) == nSamples)
      pData <- DataFrame(sampleNames)
    } else {
      stopifnot(length(sampleNames) == nSamples)
      pData <- DataFrame(sampleName=sampleNames)
    } 
    # }}}

  } else {

    # {{{ assign sampleNames  
    sampleNames <- paste0("sample", seq_len(nSamples))
    betacols <- paste0(sampleNames, ".beta")
    covgcols <- paste0(sampleNames, ".covg")
    pData <- DataFrame(sampleName=sampleNames)
    colnames(preamble) <- colNames
    # }}} 

  }
  rownames(pData) <- sampleNames
  
  nlines <- countTabix(tbx)[[1]]
  message(filename," has ",nlines," indexed loci.")
  passes <- ceiling(nlines / yieldSize)
  if (passes > 1) {
    message(filename," will require ",passes," passes of ",yieldSize," loci.")
  }
  message(filename, " looks valid for import.")

  params <- list(tbx=tbx,
                 preamble=preamble,
                 nSamples=nSamples,
                 sampleNames=sampleNames,
                 colNames=colNames,
                 betacols=betacols,
                 covgcols=covgcols,
                 hasHeader=hasHeader,
                 nlines=nlines,
                 pData=pData,
                 yieldSize=yieldSize)
  return(params) 

}
