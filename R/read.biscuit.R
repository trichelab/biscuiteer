#' a bsseq loader for Biscuit output (BED-like format, 2 or 3 cols/sample)
#' 
#' e.g. P01-028-T06.joint.ch.hg19.bed.gz has 3 samples, and 9 columns total,
#' while P01-028-T06.joint.cg.merged.hg19.bed.gz has 3 samples and 12 columns.
#' Note: the defaults assume alignment against hg19 (use genome=xyz to override)
#' Note 2: if a BED has no header, a VCF header can be used to autodetect names.
#'
#' @param BEDfile    the file (compressed or not, doesn't matter) to load
#' @param VCFfile     the file (compressed and tabixed, with header) to load
#' @param sampleNames if NULL, create; if VCF, read; if data.frame, make pData
#' @param simplify    simplify sample names by dropping .foo.bar.hg19 or similar
#' @param genome      what genome assembly were the runs aligned against? (hg19)
#' @param how         how to load the data? "data.table" (default) or "readr"
#' @param hdf5        make the object HDF5-backed? (FALSE; use in-core storage) 
#' @param sparse      are there a lot of zero-coverage sites? (TRUE, usually)
#' @param merged      are CpG sites merged? (default NULL; figure out from BED)
#' @param chunkSize   number of rows before readr reading becomes chunked (1e6)
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
read.biscuit <- function(BEDfile, 
                         VCFfile=NULL, 
                         sampleNames=NULL, 
                         simplify=FALSE, 
                         genome="hg19",
                         how=c("data.table","readr"),
                         hdf5=FALSE, 
                         sparse=TRUE,
                         merged=NULL, 
                         chunkSize=1e6, 
                         chr=NULL) { 

  how <- match.arg(how)
  params <- checkBiscuitBED(BEDfile=BEDfile, VCFfile=VCFfile, how=how, chr=chr,
                            sampleNames=sampleNames, chunk=chunkSize, hdf5=hdf5,
                            merged=merged)
  message("Reading ", ifelse(params$merged, "merged", "unmerged"), 
          " input from ", params$tbx$path, "...")

  if (params$how == "data.table") {
    # {{{
    select <- grep("\\.context", params$colNames, invert=TRUE)
    tbl <- fread(paste("gunzip -c", params$tbx$path), # for mac compatibility
                 sep="\t", sep2=",", fill=TRUE, na.string=".", select=select)
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
    # }}}
  }

  # shift from 0-based to 1-based coordinates  
  tbl[, 2] <- tbl[, 2] + 1 # FIXME: can this be done automagically? 
  message("Loaded ", params$tbx$path, ". Creating bsseq object...")
  if (params$hdf5) {
    res <- .addGenome(makeBSseq_hdf5(tbl, params, simplify=simplify), genome)
  } else { 
    res <- .addGenome(makeBSseq(tbl, params, simplify=simplify), genome)
  }
  metadata(res)$vcfHeader <- params$vcfHeader
  return(res)

}


#' @export
load.biscuit <- read.biscuit
