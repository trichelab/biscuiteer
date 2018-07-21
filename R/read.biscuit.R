#' a bsseq loader for Biscuit output (BED-like format, 2 or 3 cols/sample)
#' e.g. P01-028-T06.joint.ch.hg19.bed.gz has 3 samples, and 9 columns total,
#' while P01-028-T06.joint.cg.merged.hg19.bed.gz has 3 samples and 12 columns
#'
#' @param filename    the file (compressed or not, doesn't matter) to load
#' @param sampleNames sample names (if NULL, create; if data.frame, make pData)
#' @param hdf5        make the object HDF5-backed? (FALSE; use in-core storage) 
#' @param sparse      are there a lot of zero-coverage sites? (TRUE, usually)
#' @param chunkSize   number of rows before reading becomes chunked (1e6)
#' @param how         how to load the data? "data.table" (default) or "readr"
#' @param chr         load a specific chromosome (to rbind() later)? (NULL)
#' 
#' @return            a bsseq::BSseq object, possibly Matrix- or HDF5-backed
#'
#' @import data.table
#' @import tibble
#' @import readr
#' @import bsseq
#'
#' @seealso BSseq
#' @aliases load.biscuit
#' @seealso checkBiscuitBED
#'
#' @export
read.biscuit <- function(filename, 
                         sampleNames=NULL, 
                         hdf5=FALSE, 
                         sparse=TRUE,
                         chunkSize=1e6, 
                         how=c("data.table","readr"),
                         chr=NULL) {

  how <- match.arg(how)
  params <- checkBiscuitBED(filename,
                            sampleNames,
                            hdf5=hdf5, 
                            chunk=chunkSize, 
                            how=how,
                            chr=chr)
  message("Reading ", ifelse(params$merged, "merged", "unmerged"), 
          " input from ", params$tbx$path, "...")

  if (params$how == "readr") {
    if (params$passes > 1) { 
      f <- function(x, pos) {
        message("Reading line ", pos, "...")
        return(x)
      }
      message("Making ",params$passes," passes of ",chunkSize," loci each...")
      tbl <- with(params,
                  read_tsv_chunked(tbx$path, DataFrameCallback$new(f), na=".",
                                   skip=as.numeric(params$hasHeader), 
                                   col_names=colNames, col_types=colSpec, 
                                   chunk_size=chunkSize))
    } else { 
      message("If the following is slow, you may need to decrease chunkSize")
      message("from ",chunkSize," to something smaller & do multiple passes.")
      tbl <- with(params,
                  read_tsv(tbx$path, na=".", comment="#",
                           skip=as.numeric(params$hasHeader), 
                           col_names=colNames, col_types=colSpec))
    }
  } else if (params$how == "data.table") {
    tbl <- fread(.fixInput(params$tbx$path), sep="\t", sep2=",", fill=TRUE,
                 na.string=".", colClasses=params$colClasses)
    names(tbl) <- sub("^#", "", names(tbl))
  }

  # shift from 0-based to 1-based coordinates  
  tbl[, 2] <- tbl[, 2] + 1 
  message("Loaded ", params$tbx$path, ". Creating bsseq object...")

  if (params$hdf5) { 
    makeBSseq_hdf5(tbl, params)
  } else { 
    makeBSseq(tbl, params)
  } 

}


#' @export
load.biscuit <- read.biscuit


# helper fn
.fixInput <- function(x, mac=FALSE) { 
  if (grepl("gz$", x)) return(paste(ifelse(mac, "gzcat", "zcat"), x))
  if (grepl("bz2$", x)) return(paste("bzcat", x))
  return(x)
} 
