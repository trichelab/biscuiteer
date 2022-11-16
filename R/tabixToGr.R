#' Read from tabix-indexed bed file to list objects
#' 
#' @param paths path(s) to the bed files
#' @param chr  chromosome name
#' @param start   start coordinate of region of interest
#' @param end   end coordinate of region of interest
#' @param sample_names sample names, just use paths if not specified
#' @param is.epibed whether the input is epibed format
#' @param BPPARAM how to parallelize
#' @return a list object with DNA methylation level and depth
#'
#' @import BiocParallel
#'
#' @export
#' 

tabixRetrieve <- function(paths,
                          chr,
                          start = 1,
                          end = 2^28,
                          sample_names = NULL,
                          is.epibed = FALSE,
                          BPPARAM = SerialParam()) {
  
    input_range <- GRanges(chr, IRanges(start, end))
    df_list <- bplapply(
        paths,
        function(path) {
            df <- as.data.frame(
                t(simplify2array(strsplit(scanTabix(path, param=input_range)[[1]], '\t'))),
                stringsAsFactors = FALSE
            )
            if (is.epibed) {
                if (ncol(df) < 7 || ncol(df) > 9) stop("ERROR: Improperly formatted epiBED.")
                if (ncol(df) == 7) { # handle old format
                    colnames(df) <- c("chr", "start", "end",
                                    "readname", "read", "strand",
                                    "CG_RLE")
                    df$GC_RLE = '.'
                    df$VAR_RLE = '.'
                } else if (ncol(df) == 8) { # handle old format
                    colnames(df) <- c("chr", "start", "end",
                                    "readname", "read", "strand",
                                    "CG_RLE", "GC_RLE")
                    df$VAR_RLE = '.'
                } else { # new format
                    colnames(df) <- c("chr", "start", "end",
                                    "readname", "read", "strand",
                                    "CG_RLE", "GC_RLE", "VAR_RLE")
                }
                df[df == '.'] <- NA # replace empty strings with NA
            } else {
                colnames(df) <- c('chr','start','end','beta','depth')
            }

            df$start <- as.numeric(df$start)

            ## switch to 1-based coords
            df$start <- df$start + 1
            df$end <- as.numeric(df$end)

            # short circuit if epibed
            if (is.epibed) return(df)

            ## in case the tabix is mal-formed
            df$beta[df$beta == '.' | df$beta == 'NA'] <- NA
            df$beta <- as.numeric(df$beta)

            df$depth <- as.integer(df$depth)
            df$depth[is.na(df$beta)] <- 0

            return(df)
        }
    )

    ## make sure the coordinates are the same
    ## this is when multiple sample names or files are passed
    same_coordinates <- vapply(seq_len(length(df_list)-1),
                               function(i) identical(df_list[[i]][,seq_len(3)], df_list[[i+1]][,seq_len(3)]),
                               FUN.VALUE=logical(1))
    stopifnot(all(same_coordinates))

    ## set sample names
    if (is.null(sample_names)) {
        sample_names <- paths
    }
    stopifnot(length(sample_names) == length(paths))
    names(df_list) <- sample_names

    return(df_list)
}
