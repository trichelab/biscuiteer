#' Read in and decode the RLE representation of
#' the epibed format out of biscuit epiallele 
#' returning a read-level GRanges object.
#' 
#' NOTE: it is recommended that you ran biscuit epiallele
#' with the -g option enabled to subset to a specific region.
#' If not, this will be the equivalent of reading in the entire
#' BAM. Future work will be to subset to definied regions.
#'
#' @param epibed The path to the epibed file (must be bgzip and tabix indexed)
#' @param is.nome Whether the epibed format is derived from NOMe-seq or not
#' @param genome What genome did this come from (e.g. 'hg19')
#' @param chr Which chromosome to retrieve
#' @param start The starting position for a region of interest
#' @param end The end position for a region of interest
#'
#' @return A GRanges object
#' @export
#' 
#' @import GenomicRanges
#'
#' @examples
#' epibed.nome <- system.file("extdata", "hct116.nome.epiread.gz",
#'                            package="biscuiteer")
#' epibed.bsseq <- system.file("extdata", "hct116.bsseq.epiread.gz",
#'                             package="biscuiteer")
#' epibed.nome.gr <- readEpibed(epibed = epibed.nome, is.nome = TRUE,
#'                              genome = "hg19", chr = "chr1")
#' epibed.bsseq.gr <- readEpibed(epibed = epibed.bsseq,
#'                               genome = "hg19", chr = "chr1")

readEpibed <- function(epibed, is.nome = FALSE,
                       genome = NULL, chr = NULL,
                       start = 1, end = 2^28) {
  # check the input
  ftype <- .checkTabixFiles(epibed)
  
  if (is.null(chr)) {
    stop("Must specify chromosomes of interest.")
  }
  
  if (ftype == "multiple") {
    # read in the raw epibeds
    raw_epibed <- tabixRetrieve(epibed, chr = chr,
                                start = start,
                                end = end, is.epibed = TRUE,
                                is.nome = is.nome)
    
    # decode RLE
    # colnames should already be loaded if it's nome...
    message("Decoding RLE and converting to GRanges")
    if (is.nome) {
      epibed.gr <- lapply(raw_epibed, function(x) {
        x$CG_decode <- unlist(lapply(x$CG_RLE, .inv_rle))
        x$GC_decode <- unlist(lapply(x$GC_RLE, .inv_rle))
        x$start <- x$start + 1
        x.gr <- makeGRangesFromDataFrame(x,
                                         keep.extra.columns = TRUE)
        genome(x.gr) <- genome
        strand(x.gr) <- "*"
        return(sort(x.gr))
      })
      # make a GRangesList
      epibed.grl <- as(epibed.gr, "GRangesList")
      names(epibed.grl) <- names(raw_epibed)
    } else {
      epibed.gr <- lapply(raw_epibed, function(x) {
        x$CG_decode <- unlist(lapply(x$CG_RLE, .inv_rle))
        x$start <- x$start + 1
        x.gr <- makeGRangesFromDataFrame(x,
                                         keep.extra.columns = TRUE)
        genome(x.gr) <- genome
        strand(x.gr) <- "*"
        return(sort(x.gr))
      })
      # make a GRangesList
      epibed.grl <- as(epibed.gr, "GRangesList")
      names(epibed.grl) <- names(raw_epibed)
    }
    return(epibed.grl)
  }
  
  # in the case of a single epibed
  # read in the raw epibeds
  raw_epibed <- tabixRetrieve(epibed, chr = chr,
                              start = start,
                              end = end, is.epibed = TRUE,
                              is.nome = is.nome)
  # comes back as a list
  raw_epibed <- raw_epibed[[1]]
  message("Decoding RLE and converting to GRanges")
  # assign colnames
  if (is.nome) {
    raw_epibed$CG_decode <- unlist(lapply(raw_epibed$CG_RLE, .inv_rle))
    raw_epibed$GC_decode <- unlist(lapply(raw_epibed$GC_RLE, .inv_rle))
  } else {
    raw_epibed$CG_decode <- unlist(lapply(raw_epibed$CG_RLE, .inv_rle))
  }
  
  # shift to 1-based since this is a bed
  raw_epibed$start <- raw_epibed$start + 1
  
  # make the GRanges
  epibed_gr <- makeGRangesFromDataFrame(raw_epibed,
                                        keep.extra.columns = TRUE)
  
  # set the genome if indicated
  if (!is.null(genome)) {
    genome(epibed_gr) <- genome
  }
  
  # reset strand for now
  strand(epibed_gr) <- "*"
  
  return(sort(epibed_gr))
}

# helper
.checkTabixFiles <- function(x) {
  if (length(x) > 1) {
    # check existence 
    stopifnot(all(unlist(lapply(x, file.exists))))
    # check if gzip'd
    stopifnot(all(unlist(lapply(x, .is_gz))))
    # check for tabix indices
    stopifnot(all(unlist(lapply(x, function(x) {
      file.exists(paste0(x, ".tbi"))
    }))))
    return("multiple")
  } else {
    stopifnot(file.exists(x))
    stopifnot(.is_gz(x))
    stopifnot(file.exists(paste0(x, ".tbi")))
    return("single")
  }
}

# helper 
.inv_rle <- function(x)
{
  lengths <- as.numeric(unlist(strsplit(x, "[[:alpha:]]")))
  lengths <- lengths[-1]
  lengths[is.na(lengths)] <- 1
  lengths <- lengths[!is.na(lengths)]
  values <- unlist(strsplit(gsub("[0-9]", "", x), ""))
  values <- values[values != ""]
  uncompressed <- inverse.rle(list(lengths=lengths, values=values))
  paste(uncompressed, collapse="")
}

# helper
.is_gz <- function(f) {
  ## this is a helper from
  ## https://github.com/natverse/nat.utils/blob/master/R/gziputils.r
  if(file.access(f, mode=4)<0) return(NA)
  x=file(f)
  on.exit(close(x))
  isTRUE(summary(x)$class=='gzfile')
}

