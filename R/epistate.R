#' Tabulate to a read or fragment level matrix for
#' CGH or GCH context methylation. The input for this
#' function is the GRanges output from readEpibed().
#' It will auto-detect whether this a BS-seq or a
#' NOMe-seq run and tabulate accordingly.
#'
#' @param gr The epibed GRanges object from readEpibed()
#' 
#' @return A matrix or list of matrices
#' @export
#' 
#' @import GenomicRanges
#'
#' @examples
#' epibed.nome <- system.file("extdata", "hct116.nome.epiread",
#'                            package="biscuiteer")
#' epibed.nome.gr <- readEpibed(epibed = epibed.nome, is.nome = TRUE,
#'                              genome = "hg19")
#' epibed.tab.nome <- tabulateEpibed(epibed.nome.gr)

tabulateEpibed <- function(gr) {
  
  # check if a GRanges
  stopifnot(is(gr, "GRanges"))
  
  # autodetect if it's a NOMe epibed
  is.nome = FALSE
  if ("GC_decode" %in% names(mcols(gr))) is.nome = TRUE
  
  # tabulate strings
  # char filtered = 'F';
  # char ignored  = 'x';
  # char deletion = 'D';
  # char softclip = 'P';
  # char methylat = 'M';
  # char unmethyl = 'U';
  # char open_acc = 'O';
  # char shut_acc = 'S';
  
  # we need to build up a table of CGH and GCH
  cg_table <- .tabulateCGH(gr)
  if (is.nome) {
    gc_table <- .tabulateGCH(gr)
    return(list(cg_table = cg_table,
                gc_table = gc_table))
  }
  
  return(cg_table)
}

# helper to tabulate to a CGH table
.tabulateCGH <- function(gr) {
  # there can be duplicate read names if not collapsed to fragment
  # this can occur if reads 1 and 2 originate from the "same" strand
  # make read names unique ahead of time...
  gr$readname <- make.unique(gr$readname)
  # iterate through each read and decompose to position
  cg_gr <- do.call("c", lapply(X = 1:length(gr), FUN = function(x) {
    sub_gr <- gr[x]
    # generate a per base array
    pos_vec <- seq(start(sub_gr), end(sub_gr))
    rle_vec <- unlist(strsplit(sub_gr$CG_decode, split = ""))
    names(rle_vec) <- pos_vec
    # keep the CG status
    rle_vec_cg <- rle_vec[rle_vec %in% c("M", "U")]
    # turn back into GRanges
    rle_cg_df <- data.frame(chr = seqnames(sub_gr),
                            start = names(rle_vec_cg),
                            end = names(rle_vec_cg),
                            cg_meth_status = rle_vec_cg,
                            read_id = sub_gr$readname)
    return(makeGRangesFromDataFrame(rle_cg_df, keep.extra.columns = TRUE))
  }))
  # find the dimension of collapsed CGs to fill in the matrix
  cg_gr_len <- length(reduce(cg_gr))
  # make an empty matrix to fill in 
  cg_emp_mat <- matrix(data = NA, nrow = length(unique(cg_gr$read_id)),
                       ncol = cg_gr_len)
  rownames(cg_emp_mat) <- unique(cg_gr$read_id)
  colnames(cg_emp_mat) <- as.character(granges(reduce(cg_gr)))
  
  # go by read and extract out methylation states
  cg_gr_mat <- as.matrix(cbind(as.character(granges(cg_gr)),
                               cg_gr$read_id,
                               cg_gr$cg_meth_status))
  cg_emp_mat[cg_gr_mat[,c(2,1)]] <- cg_gr_mat[,3]
  return(cg_emp_mat)
}

# helper to tabulate to a GCH table
.tabulateGCH <- function(gr) {
  # there can be duplicate read names if not collapsed to fragment
  # this can occur if reads 1 and 2 originate from the "same" strand
  # make read names unique ahead of time...
  gr$readname <- make.unique(gr$readname)
  # iterate through each read and decompose to position
  gc_gr <- do.call("c", lapply(X = 1:length(gr), FUN = function(x) {
    sub_gr <- gr[x]
    # generate a per base array
    pos_vec <- seq(start(sub_gr), end(sub_gr))
    rle_vec <- unlist(strsplit(sub_gr$GC_decode, split = ""))
    names(rle_vec) <- pos_vec
    # keep the GC status
    rle_vec_gc <- rle_vec[rle_vec %in% c("O", "S")]
    # turn back into GRanges
    rle_gc_df <- data.frame(chr = seqnames(sub_gr),
                            start = names(rle_vec_gc),
                            end = names(rle_vec_gc),
                            gc_meth_status = rle_vec_gc,
                            read_id = sub_gr$readname)
    return(makeGRangesFromDataFrame(rle_gc_df, keep.extra.columns = TRUE))
  }))
  # find the dimension of collapsed gcs to fill in the matrix
  gc_gr_len <- length(reduce(gc_gr))
  # make an empty matrix to fill in 
  gc_emp_mat <- matrix(data = NA, nrow = length(unique(gc_gr$read_id)),
                       ncol = gc_gr_len)
  rownames(gc_emp_mat) <- unique(gc_gr$read_id)
  colnames(gc_emp_mat) <- as.character(granges(reduce(gc_gr)))
  
  # go by read and extract out methylation states
  gc_gr_mat <- as.matrix(cbind(as.character(granges(gc_gr)),
                               gc_gr$read_id,
                               gc_gr$gc_meth_status))
  gc_emp_mat[gc_gr_mat[,c(2,1)]] <- gc_gr_mat[,3]
  return(gc_emp_mat)
}

#' Plot the results of tabulateEpibed() as a 
#' quasi-lollipop plot.
#'
#' @param mat Input matrix that comes out of tabulateEpibed()
#' @param plot_read_ave Whether to also plot the average methylation state (default: TRUE)
#' @param show_readnames Whether to show the read names (default: TRUE)
#' @param show_positions Whether to show the genomic positions (default: TRUE)
#' @param meth_color What color should the methylated states be (default: 'black')
#' @param unmeth_color What color should the unmethylated states be (default: 'white')
#' @param na_color What color should the NA values be (default: 'lightgray')
#' 
#' @return An epiread ggplot object or list of ggplot objects if plot_read_ave is TRUE
#' 
#' @import ggplot2
#' 
#' @export
#'
#' @examples
#' epibed.nome <- system.file("extdata", "hct116.nome.epiread",
#'                            package="biscuiteer")
#' epibed.nome.gr <- readEpibed(epibed = epibed.nome, is.nome = TRUE,
#'                              genome = "hg19")
#' epibed.tab.nome <- tabulateEpibed(epibed.nome.gr)
#' plotEpiread(epibed.tab.nome$gc_table)

plotEpiread <- function(mat, plot_read_ave = TRUE,
                        show_readnames = TRUE,
                        show_positions = TRUE,
                        meth_color = "black",
                        unmeth_color = "white",
                        na_color = "lightgray") {
  
  # check if input is a matrix
  if (!is(mat, "matrix")) {
    stop("Input needs to be a matrix generated by tabulateEpibed()")
  }
  
  # set themes
  ql_theme <- theme_bw(12) + theme(
    axis.title = element_blank(),
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45,
                               vjust = 1,
                               hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "black"))
  
  if (!show_readnames) {
    # set the theme
    ql_theme <- ql_theme + theme(axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank())
  }
  
  if (!show_positions) {
    ql_theme <- ql_theme + theme(axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank())
  }
  
  # auto-detect input type
  is.cg = FALSE
  is.gc = FALSE
  if(any(mat %in% c("M", "U"))) is.cg = TRUE
  if(any(mat %in% c("S", "O"))) is.gc = TRUE
  
  # break if both are FALSE
  if (isFALSE(is.cg) & isFALSE(is.gc)) {
    message("Don't know what to do with input.")
    stop("Please run tabulateEpibed() first to produce input for this function.")
  }
  
  # compute average methylation state
  if (plot_read_ave) {
    # ick... copy in memory
    mat.binary <- mat
    mat.binary[mat.binary %in% c("M", "O")] <- 1
    mat.binary[mat.binary %in% c("U", "S")] <- 0
    mat.binary <- apply(mat.binary, 2, as.numeric)
    mat.meth.ave <- data.frame(ave_meth = colMeans(mat.binary, na.rm = TRUE))
    mat.meth.ave$position <- rownames(mat.meth.ave)
    mat.meth.ave$y <- "Average methylation status"
  }
  
  # cast to a 'melted' data frame
  mat.melt <- reshape2::melt(mat, id.vars = rownames(mat))
  
  # plot
  if (is.cg) {
    plt <- ggplot(mat.melt, aes(x = Var2, y = Var1)) +
      geom_point(aes(fill = value), size=6, pch=21, color="black") +
      scale_fill_manual(values = c(M=meth_color, U=unmeth_color),
                         na.value = na_color) +
      guides(color = "legend") +
      ql_theme
  } else {
    plt <- ggplot(mat.melt, aes(x = Var2, y = Var1)) +
      geom_point(aes(fill = value), size=6, pch=21, color="black") +
      scale_fill_manual(values = c(O=meth_color, S=unmeth_color),
                        na.value = na_color) +
      guides(color = "legend") +
      ql_theme
      
  }
  
  if (plot_read_ave) {
    plt_ave <- ggplot(mat.meth.ave, aes(x = position, y = y)) +
      geom_point(aes(fill = ave_meth), size=6, pch=21, color="black") +
      scale_fill_gradient(low=unmeth_color, high=meth_color) +
      guides(color = "legend") +
      ql_theme
    return(list(epistate=plt,
                meth_ave=plt_ave))
  } else {
    return(plt)
  }
}
