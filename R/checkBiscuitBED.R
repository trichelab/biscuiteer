#' inspect Biscuit BED and/or VCF output to make sure all is as it should be 
#' 
#' A BED checker for Biscuit cg/ch output (BED-like format, 2 or 3 cols/sample).
#' By default, files with over 50 million loci will be processed iteratively,
#' since data.table tends to run into problems with .gzipped joint CpH files. 
#' This function absolutely assumes that BED files are tabixed. No exceptions!
#' 
#' @param BEDfile     a BED-like file (compressed and tabixed, maybe w/header)
#' @param VCFfile     a VCF file (compressed and tabixed; only needs the header)
#' @param sampleNames if NULL create; if vector assign; if data.frame make pData
#' @param chunkSize   for files > `yieldSize` lines long, chunk the file (5e7)
#' @param merged      boolean; is this merged CpG data? (NULL; guess if merged)
#' @param hdf5        boolean; use HDF5 arrays for backing of the data? (FALSE)
#' @param sparse      boolean; use sparse Matrix objects for the data? (TRUE)
#' @param how         how to load the data? "data.table" (default) or "readr"
#' @param chr         load a specific chromosome (to rbind() later)? (NULL)
#' 
#' @return            parameters for makeBSseq or makeBSseq_hdf5
#'
#' @import            VariantAnnotation
#' @import            Rsamtools 
#' @import            readr
#'
#' @seealso           load.biscuit
#'
#' @export
checkBiscuitBED <- function(BEDfile, 
                            VCFfile=NULL, 
                            sampleNames=NULL, 
                            chunkSize=5e7, 
                            merged=NULL,
                            hdf5=FALSE,
                            sparse=TRUE,
                            how=c("data.table","readr"),
                            chr=NULL) {

  # eventual result
  params <- list() 

  # data import settings 
  params$how <- match.arg(how)
  params$sparse <- sparse
  params$hdf5 <- hdf5

  # a tabixed BED-like file that is the only mandatory argument to read.biscuit
  params$BEDfile <- BEDfile
  if (!base::grepl(".gz$", BEDfile)) stop("Only tabix'ed BEDs are supported.")
  message("Checking ", BEDfile, " for import...")
  tbx <- TabixFile(BEDfile, yieldSize=chunkSize)
  params$tbx <- tbx 

  # if we have a VCF file or VCF header, use that to get sample ordering:
  params$VCFfile <- VCFfile
  params$vcfHeader <- NULL 
  if (!is.null(params$VCFfile)) {
    if (!base::grepl(".gz$", params$VCFfile)) {
      top("Only tabix'ed VCFs are supported.")
    }
    message("Extracting sample names from ", params$VCFfile, "...") 
    params$vcfHeader <- scanVcfHeader(params$VCFfile)
    if (is.null(sampleNames)) sampleNames <- samples(params$vcfHeader)
  }

  # if restricting to one chromosome:
  params$chr <- chr 
  if (!is.null(params$chr)) {
    if (! params$chr %in% seqnamesTabix(params$tbx)) {
      stop("The requested chromosome, ", params$chr,
           ", is not found in ", params$BEDfile,".")
    } else { 
      message("Single-chromosome support is not stable yet. Beware.")
    }
  }

  # look for Tabix header, then look for problems
  params$hasHeader <- (length(headerTabix(tbx)$header) > 0)
  if (params$hasHeader) { 
    # {{{ has a header; easy and verifiable for ASM
    message(BEDfile," has a header line and is ", appendLF=FALSE)
    params$colNames <- strsplit(sub("^#","",headerTabix(tbx)$header),"\t")[[1]]
    preamble <- read.table(tbx$path, header=TRUE, sep="\t", na.strings=".", 
                           comment.char="#", nrows=3)
    colnames(preamble) <- params$colNames 
    if (is.null(merged)) merged <- any(grepl("context", colnames(preamble)))
    message(ifelse(merged, "merged", "unmerged"), " data.") 
    params$betaCols <- grep("beta", colnames(preamble), value=TRUE) 
    params$covgCols <- grep("covg", colnames(preamble), value=TRUE) 
    params$contextCols <- grep("context", colnames(preamble), value=TRUE) 
    sNames <- condenseSampleNames(tbx, stride=ifelse(merged, 3, 2))
    if (is.null(sampleNames)) sampleNames <- sNames 
    nSamples <- length(sNames)
    # }}}
  } else {
    # {{{ no header; guesswork
    message(paste0(BEDfile," has no header (!) and is "), appendLF=FALSE)
    preamble <- read.table(tbx$path, sep="\t", na.strings=".", nrows=3)
    if (is.null(merged)) merged <- grepl("(MERGE|MCG)", BEDfile, ignore=TRUE)
    message(ifelse(merged, "merged", "unmerged"), " data.") 
    cols <- c("chr","start","end")
    colsPerSample <- ifelse(merged, 3, 2)
    nSamples <- (ncol(preamble) - 3) / colsPerSample
    if (is.null(sampleNames)) sampleNames <- paste0("sample", seq_len(nSamples))
    colSuffixes <- c(".beta",".covg")
    if (merged) colSuffixes <- c(".beta",".covg",".context")
    sampcols <- paste0(rep(sampleNames, each=colsPerSample), 
                       rep(colSuffixes, nSamples))
    params$colNames <- base::gsub(" ", "", c(cols, sampcols)) # quirk
    # }}}
  }
  params$merged <- merged
  params$preamble <- preamble
  params$nSamples <- nSamples

  # by now, we know what our sampleNames are, one way or another...
  if (is(sampleNames, "DataFrame") | is(sampleNames, "data.frame")) {
    if (ncol(sampleNames) != nSamples) {
      message("The length of param$sampleNames differs from param$nSamples.")
      stop("Try setting the value of merged (TRUE or FALSE) and rerunning.")
    }
    params$pData <- DataFrame(sampleNames)
  } else {
    if (ncol(sampleNames) != nSamples) {
      message("The length of param$sampleNames differs from param$nSamples.")
      stop("Try setting the value of merged (TRUE or FALSE) and rerunning.")
    }
    params$pData <- DataFrame(sampleName=sampleNames)
  } 

  # and if they're a data.frame, use them
  if ("sampleNames" %in% names(pData)) {
    rownames(params$pData) <- pData$sampleNames  
  } else { 
    rownames(params$pData) <- params$pData[,1]
  }
  params$sampleNames <- rownames(params$pData)
  params$nLines <- countTabix(params$tbx)[[1]]
  message(BEDfile," has ", params$nLines," indexed loci.")
  if (params$how == "readr") { 
    # {{{
    params$chunkSize <- chunkSize
    params$passes <- ceiling(params$nLines / params$chunkSize)
    if(params$passes > 1) {
      message(params$BEDfile, " will require ", 
              params$passes, " passes of ", 
              params$chunkSize, " to read.")
    }
    # }}}
  }
  
  message(BEDfile, " looks valid for import.")
  params$betaCols <- grep(".beta", params$colNames, value=TRUE) 
  names(params$betaCols) <- rownames(pData)
  params$covgCols <- grep(".covg", params$colNames, value=TRUE) 
  names(params$covgCols) <- rownames(pData)
  params$contextCols <- grep(".context", params$colNames, value=TRUE) 
  names(params$contextCols) <- rownames(pData)

  # for readr 
  params$colSpec <- cols_only()
  params$colSpec[["cols"]][[params$colNames[1]]] <- col_character()
  params$colSpec[["cols"]][[params$colNames[2]]] <- col_integer()
  params$colSpec[["cols"]][[params$colNames[3]]] <- col_integer()
  for ( i in rownames(params$pData) ) { 
    params$colSpec[["cols"]][[params$betaCols[i]]] <- col_double()
    params$colSpec[["cols"]][[params$covgCols[i]]] <- col_integer()
  }
  
  # for data.table
  params$colClasses <- c() 
  params$colClasses[params$colNames[1]] <- "character"
  params$colClasses[params$colNames[2]] <- "integer"
  params$colClasses[params$colNames[3]] <- "integer"
  for ( i in rownames(params$pData) ) { 
    params$colClasses[params$betaCols[i]] <- "double"
    params$colClasses[params$covgCols[i]] <- "integer"
    params$colClasses[params$contextCols[i]]  <- "NULL"
  }
  params$colClasses <- params$colClasses[params$colNames] 
  if (params$hasHeader) {
    names(params$colClasses)[1] <- paste0("#",names(params$colClasses)[1])
  }
  return(params) 

}
