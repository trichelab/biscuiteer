#' Read in and decode the RLE representation of the epibed format out of biscuit epiread
#'
#' @param epibed The path to the epibed file (must be bgzip and tabix indexed)
#' @param genome What genome did this come from (e.g. 'hg19') (default: NULL)
#' @param chr Which chromosome to retrieve (default: NULL)
#' @param start The starting position for a region of interest (default: 1)
#' @param end The end position for a region of interest (default: 2^28)
#' @param fragment_level Whether to collapse reads to the fragment level (default: FALSE)
#'
#' @return A GRanges object
#' @export
#'
#' @import GenomicRanges
#' @import IRanges
#'
#' @examples
#'
#' epibed.nome <- system.file("extdata", "hct116.nome.epiread.gz", package="biscuiteer")
#' epibed.bsseq <- system.file("extdata", "hct116.bsseq.epiread.gz", package="biscuiteer")
#'
#' epibed.nome.gr <- readEpibed(epibed = epibed.nome, genome = "hg19", chr = "chr1")
#' epibed.bsseq.gr <- readEpibed(epibed = epibed.bsseq, genome = "hg19", chr = "chr1")
#'
readEpibed <- function(epibed,
                       genome = NULL,
                       chr = NULL,
                       start = 1,
                       end = 2^28,
                       fragment_level = FALSE) {
    # check the input
    ftype <- .checkTabixFiles(epibed)

    if (is.null(chr)) {
        stop("Must specify chromosome(s) of interest.")
    }

    if (ftype == "multiple") {
        # read in the raw epibeds
        raw_epibed <- tabixRetrieve(epibed, chr = chr, start = start, end = end, is.epibed = TRUE)

        # decode RLE
        # colnames should already be loaded...
        message("Decoding RLE and converting to GRanges")
        epibed.gr <- lapply(
            raw_epibed,
            function(x) {
                x$CG_decode <- unlist(lapply(x$CG_RLE, .inv_rle))
                x$GC_decode <- unlist(lapply(x$GC_RLE, .inv_rle))
                x$VAR_decode <- unlist(lapply(x$VAR_RLE, .inv_rle))

                x.gr <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
                genome(x.gr) <- genome
                strand(x.gr) <- "*"

                return(sort(x.gr))
            }
        )

        # make a GRangesList
        epibed.grl <- as(epibed.gr, "GRangesList")
        names(epibed.grl) <- names(raw_epibed)

        return(epibed.grl)
    }

    # in the case of a single epibed, read in the raw epibed
    raw_epibed <- tabixRetrieve(epibed, chr = chr, start = start, end = end, is.epibed = TRUE)

    # comes back as a list
    raw_epibed <- raw_epibed[[1]]
    message("Decoding RLE and converting to GRanges")

    # decode RLE strings
    raw_epibed$CG_decode  <- unlist(lapply(raw_epibed$CG_RLE, .inv_rle))
    raw_epibed$GC_decode  <- unlist(lapply(raw_epibed$GC_RLE, .inv_rle))
    raw_epibed$VAR_decode <- unlist(lapply(raw_epibed$VAR_RLE, .inv_rle))

    # make the GRanges
    epibed_gr <- makeGRangesFromDataFrame(raw_epibed, keep.extra.columns = TRUE)

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
        epibed_gr <- .collapseToFragment(epibed_gr)
    }

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
        stopifnot(all(unlist(lapply(x, function(x) { file.exists(paste0(x, ".tbi")) }))))

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
    if (is.na(x)) {
        return(NA)
    }

    lengths <- as.numeric(unlist(strsplit(x, "[[:alpha:]]")))
    lengths <- lengths[-1]
    lengths[is.na(lengths)] <- 1
    lengths <- lengths[!is.na(lengths)]
    if (!grepl("[0-9]$", x)) {
        lengths <- c(lengths, 1)
    }

    values <- unlist(strsplit(gsub("[0-9]", "", x), ""))
    values <- values[values != ""]

    uncompressed <- inverse.rle(list(lengths=lengths, values=values))
    paste(uncompressed, collapse="")
}

# helper
.recode_rle <- function(x) {
    rec_rle <- Reduce(
        paste0,
        lapply(
            seq_len(length(x$values)),
            function(l) {
                val <- x$values[l]
                len <- x$lengths[l]
                len[len == 1] <- ""
                return(paste0(val, len))
            }
        )
    )

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
.collapseToFragment <- function(gr) {
    # this will be a multi-stage collapse
    # input is a GRanges at the read level
    # assumption is that reads from the same fragment that also overlap the given region should have the same read name

    # find duplicate read names and extract
    stopifnot("readname" %in% names(mcols(gr)))
    dupe_reads <- unique(gr[duplicated(gr$readname),]$readname)
    deduped_gr <- gr[!gr$readname %in% dupe_reads,]

    # go along each read and reduce the range
    collapsed_reads <- lapply(
        dupe_reads,
        function(d) {
            dupe_pair <- gr[gr$readname == d,]
            r1 <- dupe_pair[dupe_pair$read == "1",]
            r2 <- dupe_pair[dupe_pair$read == "2",]

            ## still need to deal with padding between non-overlapping proper pairs
            if (start(r1) == start(r2)) {
                if (end(r1) >= end(r2)) {
                    # these reads perfectly overlap and should return read 1
                    r1$read <- "fragment"
                    return(r1)
                } else {
                    # this is the extreme of being considered a proper pair
                    return(.collapseCanonicalProperPair(r1, r2))
                }
            }
            if (start(r1) > start(r2)) {
                # these reads are the dovetail scenario
                # read 2 needs to be prepended to beginning of read 1
                return(.collapseDovetail(r1, r2))
            }
            if (start(r1) < start(r2)) {
                if (end(r1) >= end(r2)) {
                    r1$read <- "fragment"
                    return(r1)
                } else {
                    # these reads are prototypical
                    # read 2 needs to be appended to end of r1
                    return(.collapseCanonicalProperPair(r1, r2))
                }
            }
        }
    )

    # need to put the collapsed fragments and unpaired reads back together
    # make 'read' in mcols 'fragment'
    collapsed_reads <- do.call(c, collapsed_reads)

    # check and see if no orphan reads exist in the deduped_gr obj
    if (!length(deduped_gr)) {
      return(sort(collapsed_reads))
    }

    deduped_gr$read <- "fragment"
    collapsed_frags <- sort(c(deduped_gr, collapsed_reads))

    return(collapsed_frags)
}

# helper
.collapseDovetail <- function(r1, r2) {
    # first, check if end r2 > start r1
    # or end r2 == start r1
    cg = .mergeDovetailReads(r1, r2, "CG_decode")
    gc = .mergeDovetailReads(r1, r2, "GC_decode")
    vr = .mergeDovetailReads(r1, r2, "VAR_decode")

    # NOTE: we are in 0-based land if collapsing
    collapsed_frag <- GRanges(seqnames = seqnames(r1),
                              ranges = IRanges(start = start(r2), end = end(r1)),
                              strand = "*")

    mcols(collapsed_frag) <- data.frame(readname = r1$readname,
                                        read = "fragment",
                                        bsstrand = r1$bsstrand,
                                        CG_RLE = cg[["rle"]],
                                        GC_RLE = gc[["rle"]],
                                        VAR_RLE = vr[["rle"]],
                                        CG_decode = cg[["frag"]],
                                        GC_decode = gc[["frag"]],
                                        VAR_decode = vr[["frag"]])
    return(collapsed_frag)
}

# helper
.mergeDovetailReads <- function(r1, r2, col) {
    if (is.na(mcols(r1)[,col]) | is.na(mcols(r2)[,col])) {
        return(list("frag" = NA, "rle" = NA))
    }

    if (end(r2) >= start(r1)) {
        w <- start(r1) - start(r2)
        r2_decode <- .split_rle(mcols(r2)[,col])
        r2_decode <- r2_decode[seq_len(w)]
        r1_decode <- .split_rle(mcols(r1)[,col])
        out_frag <- paste0(Reduce(paste0, r2_decode), mcols(r1)[,col])

        # re-encode the RLE string
        out_rle <- rle(c(r2_decode, r1_decode))
    } else {
        # we are now in a situation where it's a true proper pair,
        # likely originating from same strand
        # find padding distance
        padding <- rep("x", start(r1) - (end(r2) + 1))
        r2_decode <- .split_rle(mcols(r2)[,col])
        r1_decode <- .split_rle(mcols(r1)[,col])
        out_frag <- paste0(mcols(r2)[,col], Reduce(paste0, padding), mcols(r1)[,col])

        # re-encode the RLE string
        out_rle <- rle(c(r2_decode, padding, r1_decode))
    }

    out_rle <- .recode_rle(out_rle)

    list("frag" = out_frag, "rle" = out_rle)
}

# helper
.collapseCanonicalProperPair <- function(r1, r2) {
    cg <- .mergeCanonicalProperPairReads(r1, r2, "CG_decode")
    gc <- .mergeCanonicalProperPairReads(r1, r2, "GC_decode")
    vr <- .mergeCanonicalProperPairReads(r1, r2, "VAR_decode")

    # NOTE: we are in 0-based land if collapsing
    collapsed_frag <- GRanges(seqnames = seqnames(r1),
                              ranges = IRanges(start = start(r1), end = end(r2)),
                              strand = "*")

    mcols(collapsed_frag) <- data.frame(readname = r1$readname,
                                        read = "fragment",
                                        bsstrand = r1$bsstrand,
                                        CG_RLE = cg[["rle"]],
                                        GC_RLE = gc[["rle"]],
                                        VAR_RLE = vr[["rle"]],
                                        CG_decode = cg[["frag"]],
                                        GC_decode = gc[["frag"]],
                                        VAR_decode = vr[["frag"]])

    return(collapsed_frag)
}

# helper
.mergeCanonicalProperPairReads <- function(r1, r2, col) {
    if (is.na(mcols(r1)[,col]) | is.na(mcols(r2)[,col])) {
        return(list("frag" = NA, "rle" = NA))
    }

    # the add 1 is necessary as the end pos of read 1
    # overlaps read 2 and should be filtered and thus
    # not included in the appended RLE from read 2
    w <- end(r2) - end(r1)

    if ((end(r1) - start(r2)) < 0) {
        # need to add padding
        # again, need offset to not include the end coordinate
        padding <- rep("x", start(r2) - (end(r1) + 1))
        r2_decode <- .split_rle(mcols(r2)[,col])
        r1_decode <- .split_rle(mcols(r1)[,col])
        out_frag <- paste0(mcols(r1)[,col], Reduce(paste0, padding), Reduce(paste0, r2_decode))
    } else {
        r2_decode <- .split_rle(mcols(r2)[,col])
        r2_decode <- r2_decode[(length(r2_decode) - w + 1):length(r2_decode)]
        r1_decode <- .split_rle(mcols(r1)[,col])
        out_frag <- paste0(mcols(r1)[,col], Reduce(paste0, r2_decode))
    }

    # re-encode the RLE string
    if ((end(r1) - start(r2)) < 0) {
        out_rle <- rle(c(r1_decode, padding, r2_decode))
    } else {
        out_rle <- rle(c(r1_decode, r2_decode))
    }

    out_rle <- .recode_rle(out_rle)

    list("frag" = out_frag, "rle" = out_rle)
}
