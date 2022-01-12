#' Read in and decode the RLE representation of
#' the epibed format out of biscuit epiallele 
#' returning a read-level GRanges object.
#' 
#'
#' @param epibed The path to the epibed file (must be bgzip and tabix indexed)
#' @param is.nome Whether the epibed format is derived from NOMe-seq or not
#' @param genome What genome did this come from (e.g. 'hg19')
#' @param chr Which chromosome to retrieve
#' @param start The starting position for a region of interest
#' @param end The end position for a region of interest
#' @param fragment_level Whether to collapse reads to the fragment level
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
                       start = 1, end = 2^28,
                       fragment_level = FALSE) {
  # check the input
  ftype <- .checkTabixFiles(epibed)
  
  if (is.null(chr)) {
    stop("Must specify chromosome(s) of interest.")
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
  
  # make the GRanges
  epibed_gr <- makeGRangesFromDataFrame(raw_epibed,
                                        keep.extra.columns = TRUE)
  
  # set the genome if indicated
  if (!is.null(genome)) {
    genome(epibed_gr) <- genome
  }
  
  # reset strand for now
  strand(epibed_gr) <- "*"
  
  # collapse to fragment
  if (fragment_level) {
    message("Collapsing to fragment level")
    message("This will take some time if a large region is being analyzed")
    epibed_gr <- .collapseToFragment(epibed_gr,
                                     is.nome = is.nome)
  }
  
  # shift to 1-based since this is a bed
  start(epibed_gr) <- start(epibed_gr) + 1
  
  # make sure it's sorted
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
.inv_rle <- function(x) {
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
.recode_rle <- function(x) {
  rec_rle <- Reduce(paste0, lapply(1:length(x$values), function(l) {
    val <- x$values[l]
    len <- x$lengths[l]
    len[len == 1] <- ""
    return(paste0(val, len))
  }))
  return(rec_rle)
}

# helper
.split_rle <- function(x) {
  return(unlist(strsplit(x, split = "")))
}

# helper
.is_gz <- function(f) {
  # this is a helper from
  # https://github.com/natverse/nat.utils/blob/master/R/gziputils.r
  if(file.access(f, mode=4)<0) return(NA)
  x=file(f)
  on.exit(close(x))
  isTRUE(summary(x)$class=='gzfile')
}

# helper
.collapseToFragment <- function(gr, is.nome = FALSE) {
  # input is a GRanges at the read level
  # assumption is that reads from the same fragment
  # that also overlap the given region
  # should have the same read name
  # this will be a multi-stage collapse
  # find duplicate read names and extract
  stopifnot("readname" %in% names(mcols(gr)))
  dupe_reads <- unique(gr[duplicated(gr$readname),]$readname)
  deduped_gr <- gr[!gr$readname %in% dupe_reads,]
  # go along each read and reduce the range
  collapsed_reads <- lapply(dupe_reads, function(d) {
    dupe_pair <- gr[gr$readname == d,]
    r1 <- dupe_pair[dupe_pair$read == "1",]
    r2 <- dupe_pair[dupe_pair$read == "2",]
    ## still need to deal with padding between non-overlapping proper pairs
    if (start(r1) == start(r2) & end(r1) == end(r2)) {
      # these reads perfectly overlap and should return read 1
      r1$read <- "fragment"
      return(r1)
    }
    if (start(r1) > start(r2)) {
      # these reads are the dovetail scenario
      # read 2 needs to be appended to beginning of read 1
      return(.collapseDovetail(r1, r2,
                               is.nome = is.nome))
    }
    if (start(r1) < start(r2)) {
      # these reads are prototypical
      # read 2 needs to be appended to end of r1
      return(.collapseProperPair(r1, r2,
                                 is.nome = is.nome))
    }
  })
  ## need to put the collapsed fragments and unpaired reads
  ## back together
  ## make 'read' in mcols 'fragment'
  collapsed_reads <- do.call(c, collapsed_reads)
  deduped_gr$read <- "fragment"
  collapsed_frags <- sort(c(deduped_gr,
                            collapsed_reads))
  return(collapsed_frags)
}

# helper
.collapseDovetail <- function(r1, r2,
                              is.nome = FALSE) {
  # first, check if end r2 > start r1
  if (end(r2) > start(r1)) {
    w <- start(r1) - start(r2)
    r2_decode <- .split_rle(r2$CG_decode)
    r2_decode <- r2_decode[1:w]
    r1_decode <- .split_rle(r1$CG_decode)
    cg_frag <- paste0(Reduce(paste0, r2_decode),
                      r1$CG_decode)
    # re-encode the RLE string
    cg_rle <- rle(c(r2_decode, r1_decode))
  } else {
    # we are now in a situation where it's a true
    # proper pair, likely originating from same strand
    # find padding distance
    padding <- rep("x", start(r1) - (end(r2) + 1))
    r2_decode <- .split_rle(r2$CG_decode)
    r1_decode <- .split_rle(r1$CG_decode)
    cg_frag <- paste0(r2$CG_decode,
                      Reduce(paste0, padding),
                      r1$CG_decode)
    # re-encode the RLE string
    cg_rle <- rle(c(r2_decode,
                    padding,
                    r1_decode))
  }
  
  cg_rle <- .recode_rle(cg_rle)
  # NOTE: we are in 0-based land if collapsing
  collapsed_frag <- GRanges(seqnames = seqnames(r1),
                            ranges = IRanges(start = start(r2),
                                             end = end(r1)),
                            strand = "*")
  if (is.nome) {
    if (end(r2) > start(r1)) {
      r2_decode <- .split_rle(r2$GC_decode)
      r2_decode <- r2_decode[1:w]
      r1_decode <- .split_rle(r1$GC_decode)
      gc_frag <- paste0(Reduce(paste0, r2_decode),
                        r1$GC_decode)
      # re-encode the RLE string
      gc_rle <- rle(c(r2_decode, r1_decode))
    } else {
      # we are now in a situation where it's a true
      # proper pair, likely originating from same strand
      r2_decode <- .split_rle(r2$GC_decode)
      r1_decode <- .split_rle(r1$GC_decode)
      gc_frag <- paste0(r2$GC_decode,
                        Reduce(paste0, padding),
                        r1$GC_decode)
      # re-encode the RLE string
      gc_rle <- rle(c(r2_decode,
                      padding,
                      r1_decode))
    }
    gc_rle <- .recode_rle(gc_rle)
    collapsed_frag.meta <- data.frame(readname = r1$readname,
                                      read = "fragment",
                                      CG_RLE = cg_rle,
                                      GC_RLE = gc_rle,
                                      CG_decode = cg_frag,
                                      GC_decode = gc_frag)
  } else {
    collapsed_frag.meta <- data.frame(readname = r1$readname,
                                      read = "fragment",
                                      CG_RLE = cg_rle,
                                      CG_decode = cg_frag)
  }
  mcols(collapsed_frag) <- collapsed_frag.meta
  return(collapsed_frag)
}

# helper
.collapseProperPair <- function(r1, r2,
                                is.nome = is.nome) {
  # the add 1 is necessary as the end pos of read 1
  # overlaps read 2 and should be filtered and thus
  # not included in the appended RLE from read 2
  w <- end(r2) - end(r1)
  if ((end(r1) - start(r2)) < 0) {
    # need to add padding
    # again, need to offset to not include the end coordinate
    padding <- rep("x", start(r2) - (end(r1) + 1))
    r2_decode <- .split_rle(r2$CG_decode)
    r1_decode <- .split_rle(r1$CG_decode)
    cg_frag <- paste0(r1$CG_decode,
                      Reduce(paste0, padding),
                      Reduce(paste0, r2_decode))
  } else {
    r2_decode <- .split_rle(r2$CG_decode)
    r2_decode <- r2_decode[(length(r2_decode) - w + 1):length(r2_decode)]
    r1_decode <- .split_rle(r1$CG_decode)
    cg_frag <- paste0(r1$CG_decode,
                      Reduce(paste0, r2_decode))
  }
  
  # re-encode the RLE string
  if ((end(r1) - start(r2)) < 0) {
    cg_rle <- rle(c(r1_decode,
                    padding,
                    r2_decode))
  } else {
    cg_rle <- rle(c(r1_decode, r2_decode))
  }
  cg_rle <- .recode_rle(cg_rle)
  # NOTE: we are in 0-based land if collapsing
  collapsed_frag <- GRanges(seqnames = seqnames(r1),
                            ranges = IRanges(start = start(r1),
                                             end = end(r2)),
                            strand = "*")
  if (is.nome) {
    if ((end(r1) - start(r2)) < 0) {
      r2_decode <- .split_rle(r2$GC_decode)
      r1_decode <- .split_rle(r1$GC_decode)
      gc_frag <- paste0(r1$GC_decode,
                        Reduce(paste0, padding),
                        Reduce(paste0, r2_decode))
      gc_rle <- rle(c(r1_decode,
                      padding,
                      r2_decode))
    } else {
      r2_decode <- .split_rle(r2$GC_decode)
      r2_decode <- r2_decode[(length(r2_decode) - w + 1):length(r2_decode)]
      r1_decode <- .split_rle(r1$GC_decode)
      gc_frag <- paste0(r1$GC_decode,
                        Reduce(paste0, r2_decode))
      gc_rle <- rle(c(r1_decode, r2_decode))
    }
    gc_rle <- .recode_rle(gc_rle)
    collapsed_frag.meta <- data.frame(readname = r1$readname,
                                      read = "fragment",
                                      CG_RLE = cg_rle,
                                      GC_RLE = gc_rle,
                                      CG_decode = cg_frag,
                                      GC_decode = gc_frag)
  } else {
    collapsed_frag.meta <- data.frame(readname = r1$readname,
                                      read = "fragment",
                                      CG_RLE = cg_rle,
                                      CG_decode = cg_frag)
  }
  mcols(collapsed_frag) <- collapsed_frag.meta
  return(collapsed_frag)
}

#' Function to remove insertions to reset genomic coordinates
#' for plotting of read or fragment level methylation states,
#' relative to the reference genome.
#'
#' @param gr Input is a read or fragment level GRanges out of readEpibed
#' @param is.nome Whether the input is NOMe-seq or not
#'
#' @return A new GRanges object with reset genomic coordinates and filtered
#' CG and/or GC RLE strings
#' 
#' @import GenomicRanges
#' @export
#'
#' @examples
#' 
#' epibed.nome <- system.file("extdata", "hct116.nome.epiread.gz",
#'                            package="biscuiteer")
#' epibed.bsseq <- system.file("extdata", "hct116.bsseq.epiread.gz",
#'                             package="biscuiteer")
#' epibed.nome.gr <- readEpibed(epibed = epibed.nome, is.nome = TRUE,
#'                              genome = "hg19", chr = "chr1")
#' epibed.bsseq.gr <- readEpibed(epibed = epibed.bsseq,
#'                               genome = "hg19", chr = "chr1")
#'                               
#' epibed.nome.filt.gr <- removeInsertions(epibed.nome.gr, is.nome = TRUE)
#' epibed.bsseq.filt.gr <- removeInsertions(epibed.bsseq.gr)

removeInsertions <- function(gr, is.nome = FALSE) {
  # wrapper to remove insertions from reads
  filt_gr <- do.call(c, lapply(1:length(gr), function(x) {
    sub_gr <- gr[x]
    # NOTE: input is 1-based
    start(sub_gr) <- start(sub_gr) - 1
    # generate a per base array
    pos_vec <- seq(start(sub_gr), end(sub_gr))
    if (is.nome) {
      rle_vec.cg <- .split_rle(sub_gr$CG_decode)
      rle_vec.gc <- .split_rle(sub_gr$GC_decode)
    } else {
      rle_vec.cg <- .split_rle(sub_gr$CG_decode)
    }
    names(rle_vec.cg) <- pos_vec
    if (is.nome) names(rle_vec.gc) <- pos_vec
    # filter out insertions
    rle_vec.cg <- .filterInsertions(rle_vec.cg)
    if (is.nome) rle_vec.gc <- .filterInsertions(rle_vec.gc)
    # check
    if (is.nome) stopifnot(all(names(rle_vec.cg) == names(rle_vec.gc)))
    # add back the new, filtered fragments
    sub_gr$CG_decode <- Reduce(paste0, rle_vec.cg)
    if (is.nome) sub_gr$GC_decode <- Reduce(paste0, rle_vec.gc)
    # reset the RLE string and genomic coords
    # back to 1-based
    start(sub_gr) <- as.integer(names(rle_vec.cg)[1]) + 1
    end(sub_gr) <- as.integer(names(rle_vec.cg)[length(rle_vec.cg)])
    sub_gr$CG_RLE <- .recode_rle(rle(.split_rle(sub_gr$CG_decode)))
    if (is.nome) sub_gr$GC_RLE <- .recode_rle(rle(.split_rle(sub_gr$GC_decode)))
    return(sub_gr)
  }))
  return(filt_gr)
}

# helper to filter out indels
# this is needed to reset coordinates properly
.filterInsertions <- function(readlvl_vec) {
  # the input here is a named vec of positions
  # we need to pull everything out that is not
  # a lower case a,c,g,t
  # we can correct for new starts if someone has
  # not filtered the first few bases
  exclude_bases <- c("a", "c",
                     "g", "t")
  # note: if a SNP has a base, it will be upper case
  # case sensitivity matters here
  # grab the start from the original string
  strt <- names(readlvl_vec)[1]
  filtrd_vec <- readlvl_vec[!readlvl_vec %in% exclude_bases]
  # short circuit if nothing is filtered
  if (suppressWarnings(all(names(filtrd_vec) == names(readlvl_vec)))) {
    return(readlvl_vec)
  }
  # if the first base is no longer equal to the start
  # from original read, reset to new start
  if (strt != names(filtrd_vec)[1]) {
    strt <- names(filtrd_vec)[1]
  }
  names(filtrd_vec) <- seq(as.numeric(strt),
                           c(as.numeric(strt)+length(filtrd_vec)-1))
  return(filtrd_vec)
}
