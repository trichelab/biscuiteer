#' Read biscuit output into bsseq object
#'
#' Takes BED-like format with 2 or 3 columns per sample. Unmerged CpG files
#' have 2 columns (beta values and coverage), whereas merged CpG files have
#' 3 columns (beta values, coverage, and context).
#'
#' NOTE: Assumes alignment against hg19 (use genome argument to override).
#' NOTE: Requires header from VCF file to detect sample names
#'
#' @param BEDfile      A BED-like file - must be compressed and tabix'ed
#' @param VCFfile      A VCF file - must be compressed and tabix'ed. Only the
#'                       header information is needed.
#' @param merged       Is this merged CpG data?
#' @param sampleNames  Names of samples - NULL: create names, vector: assign
#'                       names, data.frame: make pData (DEFAULT: NULL)
#' @param simplify     Simplify sample names by dropping .foo.bar.hg19? (or
#'                       similar) (DEFAULT: FALSE)
#' @param genome       Genome assembly the runs were aligned against
#'                       (DEFAULT: "hg19")
#' @param how          How to load data - either data.table or readr
#'                       (DEFAULT: "data.table")
#' @param hdf5         Make the object HDF5-backed - CURRENTLY NOT AVAILABLE
#'                       (DEFAULT: FALSE)
#' @param hdf5dir      Directory to store HDF5 files if 'hdf5' = TRUE
#'                       (DEFAULT: NULL)
#' @param sparse       Use sparse Matrix objects for the data? (DEFAULT: FALSE)
#' @param chunkSize    Number of rows before readr reading becomes chunked
#'                       (DEFAULT: 1e6)
#' @param chr          Load a specific chromosome? (DEFAULT: NULL)
#' @param which        A GRanges of regions to load - NULL loads them all
#'                       (DEFAULT: NULL)
#' @param verbose      Print extra statements? (DEFAULT: FALSE)
#'
#' @return             A bsseq::BSseq object
#'
#' @importFrom data.table fread
#' @importFrom R.utils gunzip
#' @importFrom rtracklayer export
#' @import SummarizedExperiment
#' @import readr
#' @import bsseq
#'
#' @seealso bsseq
#' @seealso checkBiscuitBED
#'
#' @aliases loadBiscuit
#'
#' @examples
#'
#'   orig_bed <- system.file("extdata", "MCF7_Cunha_chr11p15.bed.gz",
#'                           package="biscuiteer")
#'   orig_vcf <- system.file("extdata", "MCF7_Cunha_header_only.vcf.gz",
#'                           package="biscuiteer")
#'   bisc <- readBiscuit(BEDfile = orig_bed, VCFfile = orig_vcf,
#'                       merged = FALSE)
#'
#' @export
#'
readBiscuit <- function(BEDfile, 
                        VCFfile, 
                        merged, 
                        sampleNames = NULL, 
                        simplify = FALSE, 
                        genome = "hg19",
                        how = c("data.table", "readr"),
                        hdf5 = FALSE, 
                        hdf5dir = NULL,
                        sparse = FALSE,
                        chunkSize = 1e6, 
                        chr = NULL,
                        which = NULL,
                        verbose = FALSE) { 

  # Check if required inputs are missing
  # Print more useful messages if they are
  if (rlang::is_missing(BEDfile))
    stop("Tabix'ed BED file from biscuit is required.\n")
  if (rlang::is_missing(VCFfile)) {
    stop("Tabix'ed VCF file from biscuit is required. Header information is ",
         "used to set up column names.")
  }
  if (rlang::is_missing(merged)) {
    stop("merged flag is required. merged = TRUE if 'biscuit mergecg' was run ",
         "after 'biscuit vcf2bed'. Otherwise use merged = FALSE.")
  }

  how <- match.arg(how)
  params <- checkBiscuitBED(BEDfile=BEDfile, VCFfile=VCFfile, how=how, chr=chr,
                            sampleNames=sampleNames, chunkSize=chunkSize,
                            hdf5=hdf5, sparse=sparse, merged=merged)
  message("Reading ", ifelse(params$merged, "merged", "unmerged"), 
          " input from ", params$tbx$path, "...")

  if (params$how == "data.table") {
    # {{{
    select <- grep("\\.context", params$colNames, invert=TRUE)
    if (is.null(which)) {
      tbl <- fread(gunzip(params$tbx$path, remove = FALSE), sep="\t", sep2=",",
                   fill=TRUE, na.strings=".", select=select)
      unzippedName <- sub("\\.gz$", "", params$tbx$path)
      if (file.exists(unzippedName)) {
        file.remove(unzippedName)
      }
    } else { 
      tmpBed <- tempfile(fileext=".bed")
      export(which, tmpBed)
      cmd <- paste("tabix -R", tmpBed, params$tbx$path)
      tbl <- fread(cmd=cmd, sep="\t", sep2=",",
                   fill=TRUE, na.strings=".", select=select)
    }
    if (params$hasHeader == FALSE) names(tbl) <- params$colNames[select]
    names(tbl) <- sub("^#", "", names(tbl))
    # }}}
  } else if (params$how == "readr") {
    # {{{
    if (params$passes > 1) { 
      f <- function(x, pos) {
        message("Reading line ", pos, "...")
        return(x)
      }
      message("Making ",params$passes," passes of ",chunkSize," loci each...")
      tbl <- read_tsv_chunked(params$tbx$path, DataFrameCallback$new(f), na=".",
                              skip=as.numeric(params$hasHeader), 
                              col_names=params$colNames,
                              col_types=params$colSpec, chunk_size=chunkSize)
    } else { 
      message("If the following is slow, you may need to decrease chunkSize")
      message("from ",chunkSize," to something smaller & do multiple passes.")
      tbl <- read_tsv(params$tbx$path, na=".", comment="#",
                      skip=as.numeric(params$hasHeader), 
                      col_names=params$colNames, col_types=params$colSpec)
    }
    # }}}
  }

  # shift from 0-based to 1-based coordinates  
  tbl[, 2] <- tbl[, 2] + 1 # FIXME: can this be done automagically?

  # Remove CpG sites with zero-coverage
  if(!params$sparse) {
    message("Excluding CpG sites with uniformly zero coverage...")
    tbl <- tbl[rowSums(is.na(tbl)) == 0, ]
  }

  message("Loaded ", params$tbx$path, ". Creating bsseq object...")
  res <- makeBSseq(tbl, params, simplify=simplify, verbose=verbose)
  metadata(res)$vcfHeader <- params$vcfHeader
  genome(rowRanges(res)) <- genome

  if (hdf5) {
    if (is.null(hdf5dir)) {
      stop("You must provide an `hdf5dir` argument if you set `hdf5` to TRUE.")
    } else {
      if (dir.exists(hdf5dir)) {
        stop("The directory you specified already exists!")
      } else { 
        res <- HDF5Array::saveHDF5SummarizedExperiment(res, dir=hdf5dir)
      }
    }
  } 
  return(res)

}


#' @describeIn readBiscuit Alias for readBiscuit
#'
loadBiscuit <- readBiscuit
