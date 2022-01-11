#' Read in and decode the RLE representation of
#' the epibed format out of biscuit epiallele 
#' returning a read-level GRanges object.
#' 
#' NOTE: it is recommended that you ran biscuit epiallele
#' with the -g option enabled to subset to a specific region.
#' If not, this will be the equivalent of reading in the entire
#' BAM. Future work will be to subset to definied regions.
#'
#' @param epibed The path to the epibed file
#' @param is.nome Whether the epibed format is derived from NOMe-seq or not
#' @param genome What genome did this come from (e.g. 'hg19')
#'
#' @return A GRanges object
#' @export
#' 
#' @import GenomicRanges
#' @importFrom data.table fread
#' @importFrom R.utils gunzip
#'
#' @examples
#' epibed.nome <- system.file("extdata", "hct116.nome.epiread",
#'                            package="biscuiteer")
#' epibed.bsseq <- system.file("extdata", "hct116.bsseq.epiread",
#'                             package="biscuiteer")
#' epibed.nome.gr <- readEpibed(epibed = epibed.nome, is.nome = TRUE,
#'                              genome = "hg19")
#' epibed.bsseq.gr <- readEpibed(epibed = epibed.bsseq,
#'                               genome = "hg19")

readEpibed <- function(epibed, is.nome = FALSE,
                       genome = NULL) {
  # check if the file exists
  stopifnot(file.exists(epibed))
  # read in the entire file for now because we assume
  # that the user has done something sane like used
  # the -g option in biscuit epiread
  message("Reading in ", epibed)
  
  # check if gzipped
  if (.is_gz(epibed)) {
    raw_epibed <- fread(gunzip(epibed, remove = FALSE), sep="\t")
  } else {
    raw_epibed <- fread(epibed, sep="\t")
  }
  
  message("Decoding RLE and converting to GRanges")
  # assign colnames
  if (is.nome) {
    colnames(raw_epibed) <- c("chr", "start", "end",
                              "readname", "read", "strand",
                              "CG_RLE", "GC_RLE")
    raw_epibed$CG_decode <- unlist(lapply(raw_epibed$CG_RLE, .inv_rle))
    raw_epibed$GC_decode <- unlist(lapply(raw_epibed$GC_RLE, .inv_rle))
  } else {
    colnames(raw_epibed) <- c("chr", "start", "end",
                              "readname", "read", "strand",
                              "CG_RLE")
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
