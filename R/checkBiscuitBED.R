#' A BED checker for Biscuit cg/ch output (BED-like format, 2 or 3 cols/sample).
#' By default, files with over 50 million loci will be processed iteratively,
#' since data.table tends to run into problems with .gzipped joint CpH files. 
#' This function absolutely assumes that BED files are tabixed. No exceptions!
#' 
#' @param filename    the file (compressed and tabixed, with header) to load
#' @param sampleNames sample names (if NULL, create; if data.frame, make pData)
#' @param chunkSize   for files > `yieldSize` lines long, chunk the file (5e7)
#' @param merged      boolean; is this merged CpG data? (NULL; guess if merged)
#' @param hdf5        boolean; use HDF5 arrays for backing of the data? (FALSE)
#' @param sparse      boolean; use sparse Matrix objects for the data? (TRUE)
#' @param how         how to load the data? "readr" (default) or "data.table"
#' 
#' @return            parameters for makeBSseq or makeBSseq_hdf5
#'
#' @import            Rsamtools 
#' @import            readr
#'
#' @seealso           load.biscuit
#'
#' @export
checkBiscuitBED <- function(filename, 
                            sampleNames=NULL, 
                            chunkSize=5e7, 
                            merged=NULL,
                            hdf5=FALSE,
                            sparse=TRUE,
                            how=c("data.table","readr")){

  how <- match.arg(how)
  if (!base::grepl(".gz$", filename)) stop("Only tabix'ed files are supported.")
  message("Checking ", filename, " for import...")
  tbx <- TabixFile(filename, yieldSize=chunkSize)

  # look for Tabix header, then look for problems
  hasHeader <- (length(headerTabix(tbx)$header) > 0)
  if (hasHeader) { 
    # {{{ has a header; easy and verifiable for ASM
    message(filename," has a header line and is ", appendLF=FALSE)
    colNames <- strsplit(sub("^#", "", headerTabix(tbx)$header), "\t")[[1]]
    preamble <- read.table(tbx$path, header=TRUE, sep="\t", na.strings=".", 
                           comment.char="#", nrows=3)
    colnames(preamble) <- colNames 
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
    if (is.null(sampleNames)) {
      sampleNames <- paste0("sample", seq_len(nSamples))
    }
    colSuffixes <- c(".beta",".covg")
    if (merged) colSuffixes <- c(".beta",".covg",".context")
    sampcols <- paste0(rep(sampleNames, each=colsPerSample), 
                       rep(colSuffixes, nSamples))
    colNames <- base::gsub(" ", "", c(cols, sampcols)) # quirk
    # }}}
  }

  # by now, we know what our sampleNames are, one way or another...
  if (is(sampleNames, "DataFrame") | is(sampleNames, "data.frame")) {
    stopifnot(ncol(sampleNames) == nSamples)
    pData <- DataFrame(sampleNames)
  } else {
    stopifnot(length(sampleNames) == nSamples)
    pData <- DataFrame(sampleName=sampleNames)
  } 

  # and if they're a data.frame, use them
  if ("sampleNames" %in% names(pData)) {
    rownames(pData) <- pData$sampleNames  
  } else { 
    rownames(pData) <- pData[,1]
  }
  
  nlines <- countTabix(tbx)[[1]]
  message(filename," has ",nlines," indexed loci.")
  passes <- ceiling(nlines / chunkSize)
  if (passes > 1) message(filename," takes ",passes," passes of ",chunkSize,".")
  message(filename, " looks valid for import.")
  betacols <- grep(".beta", colNames, value=TRUE) 
  names(betacols) <- rownames(pData)
  covgcols <- grep(".covg", colNames, value=TRUE) 
  names(covgcols) <- rownames(pData)

  # for readr (or data.table dropcols)
  colSpec <- cols_only()
  colSpec[["cols"]][[colNames[1]]] <- col_character()
  colSpec[["cols"]][[colNames[2]]] <- col_integer()
  colSpec[["cols"]][[colNames[3]]] <- col_integer()
  for ( i in rownames(pData) ) { 
    colSpec[["cols"]][[betacols[i]]] <- col_double()
    colSpec[["cols"]][[covgcols[i]]] <- col_integer()
  }

  params <- list(tbx=tbx,
                 merged=merged,
                 preamble=preamble,
                 nSamples=nSamples,
                 sampleNames=rownames(pData),
                 colNames=colNames,
                 colSpec=colSpec,
                 betacols=betacols,
                 covgcols=covgcols,
                 hasHeader=hasHeader,
                 nlines=nlines,
                 pData=pData,
                 passes=passes,
                 chunkSize=chunkSize,
                 sparse=sparse,
                 hdf5=hdf5,
                 how=how)
  return(params) 

}
