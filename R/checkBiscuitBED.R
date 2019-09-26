#' Inspect Biscuit VCF and BED files
#'
#' A BED checker for Biscuit CpG/CpH output (BED-like format with 2 or 3
#' columns per sample). By default, files with more than 50 million loci will
#' be processed iteratively, since data.table tends to run into problems with
#' gzipped joint CpH files.
#'
#' Input BED and VCF files must be tabix'ed. No exceptions!
#'
#' @param BEDfile      A BED-like file - must be compressed and tabix'ed
#' @param VCFfile      A VCF file - must be compressed and tabix'ed. Only the
#'                       header information is needed.
#' @param merged       Is this merged CpG data?
#' @param sampleNames  Names of samples - NULL: create names, vector: assign
#'                       names, data.frame: make pData (DEFAULT: NULL)
#' @param chunkSize    For files longer than `yieldSize` number of lines long,
#'                       chunk the file (DEFAULT: 5e7)
#' @param hdf5         Use HDF5 arrays for backing the data? (DEFAULT: FALSE)
#' @param sparse       Use sparse Matrix objects for the data? (DEFAULT: TRUE)
#' @param how          How to load the data - "data.table" or "readr"?
#'                       (DEFAULT: data.table)
#' @param chr          Load a specific chromosome to rbind() later?
#'                       (DEFAULT: NULL)
#'
#' @return             Parameters to be supplied to makeBSseq
#'
#' @importFrom methods as is
#' @importFrom utils read.table
#' @import VariantAnnotation
#' @import Rsamtools
#' @import readr
#'
#' @seealso read.biscuit
#'
#' @examples
#'
#'   orig_bed <- system.file("extdata", "MCF7_Cunha_chr11p15.bed.gz",
#'                           package="biscuiteer")
#'   orig_vcf <- system.file("extdata", "MCF7_Cunha_header_only.vcf.gz",
#'                           package="biscuiteer")
#'   params <- checkBiscuitBED(BEDfile = orig_bed, VCFfile = orig_vcf,
#'                             merged = FALSE)
#'
#' @export
#'
checkBiscuitBED <- function(BEDfile, 
                            VCFfile, 
                            merged,
                            sampleNames = NULL, 
                            chunkSize = 5e7, 
                            hdf5 = FALSE,
                            sparse = TRUE,
                            how = c("data.table","readr"),
                            chr = NULL) {

  # Check if required inputs are missing
  # Print more useful messages if they are
  if (rlang::is_missing(BEDfile))
    stop("Tabix'ed BED file from biscuit is required.\n")
  if (rlang::is_missing(VCFfile)) {
    err_message <- paste("Tabix'ed VCF file from biscuit is required.",
                         "Header information is used to set up column names.\n")
    stop(err_message)
  }
  if (rlang::is_missing(merged)) {
    err_message <- paste("merged flag is required.",
                         "merged = TRUE if 'biscuit mergecg' was",
                         "run after 'biscuit vcf2bed'.",
                         "Otherwise use merged = FALSE.\n")
    stop(err_message)
  }

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
      stop("Only tabix'ed VCFs are supported.")
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

  # Set up columns based on number of samples in header of BED file or VCF file
  # Do some checks to make sure things are kosher
  params$hasHeader <- (length(headerTabix(tbx)$header) > 0)
  if (params$hasHeader) {
    # BED file has a header line ---> Not standard on a biscuit BED file
    # TODO: Can probably remove this as a biscuit BED
    #       file should never have this information
    message(BEDfile, " has a header line")

    params$colNames <- strsplit(sub("^#","",headerTabix(tbx)$header),"\t")[[1]]
    preamble <- read.table(tbx$path, header=TRUE, sep="\t", na.strings=".", 
                           comment.char="#", nrows=3)
    colnames(preamble) <- params$colNames
    
    message("Assuming ", ifelse(merged, "merged", "unmerged"),
            " data. Checking now...", appendLF=FALSE)
    if (!merged & any(grepl("context", colnames(preamble)))) {
      stop("merged flag is FALSE, but appears to be a merged BED file. ",
           "Set merged = TRUE.")
    }
    if (merged & !any(grepl("context", colnames(preamble)))) {
      stop("merged flag is TRUE, but does not appear to be a merged BED file. ",
           "Set merged = FALSE.")
    }
    message(" ...", ifelse(merged, "merged", "unmerged"), " file seems okay.")

    params$betaCols <- grep("beta", colnames(preamble), value=TRUE) 
    params$covgCols <- grep("covg", colnames(preamble), value=TRUE) 
    if (merged)
      params$contextCols <- grep("context", colnames(preamble), value=TRUE) 
    sNames <- condenseSampleNames(tbx, stride=ifelse(merged, 3, 2))
    if (is.null(sampleNames)) sampleNames <- sNames 
    nSamples <- length(sNames)
  } else {
    message(paste0(BEDfile, " does not have a header. ",
                   "Using VCF file header information ",
                   "to help set column names."))

    preamble <- read.table(tbx$path, sep="\t", na.strings=".", nrows=3)
    cols <- c("chr","start","end")
    colsPerSample <- ifelse(merged, 3, 2)
    nSamples <- (ncol(preamble) - 3) / colsPerSample

    message("Assuming ", ifelse(merged, "merged", "unmerged"),
            " data. Checking now...", appendLF=FALSE)
    if (((ncol(preamble) - 3) %% colsPerSample) != 0) {
      stop(paste0("Number of columns per sample does not seem to line up. ",
                  "Double check your merge flag is correct."))
    }
    if (nSamples != length(sampleNames)) {
        stop(paste0("nSamples (", nSamples, ") does not match the ",
                    "number of sampleNames (", length(sampleNames),")"))
    }
    message(paste0(" ...The file might be alright. ",
                   "Double check if you're worried."))

    if (is.null(sampleNames)) sampleNames <- paste0("sample", seq_len(nSamples))
    colSuffixes <- c(".beta",".covg")
    if (merged) colSuffixes <- c(".beta",".covg",".context")
    sampcols <- paste0(rep(sampleNames, each=colsPerSample), 
                       rep(colSuffixes, nSamples))
    params$colNames <- base::gsub(" ", "", c(cols, sampcols)) # quirk
  }
  params$merged <- merged
  params$preamble <- preamble
  params$nSamples <- nSamples

  # by now, we know what our sampleNames are, one way or another...
  if (is(sampleNames, "DataFrame") | is(sampleNames, "data.frame")) {
    if (nrow(sampleNames) != nSamples) {
      message("The length of param$sampleNames differs from param$nSamples.")
      stop("Try setting the value of merged (TRUE or FALSE) and rerunning.")
    }
    params$pData <- DataFrame(sampleNames)
  } else {
    if (length(sampleNames) != nSamples) {
      message("The length of param$sampleNames differs from param$nSamples.")
      stop("Try setting the value of merged (TRUE or FALSE) and rerunning.")
    }
    params$pData <- DataFrame(sampleNames=sampleNames)
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
  for ( i in seq_along(rownames(params$pData)) ) {
    params$colSpec[["cols"]][[params$betaCols[i]]] <- col_double()
    params$colSpec[["cols"]][[params$covgCols[i]]] <- col_integer()
  }
  
  # for data.table
  params$colClasses <- c() 
  params$colClasses[params$colNames[1]] <- "character"
  params$colClasses[params$colNames[2]] <- "integer"
  params$colClasses[params$colNames[3]] <- "integer"
  for ( i in seq_along(rownames(params$pData)) ) {
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
