#' Choose loci or features by extremality
#'
#' This function finds the k most extremal features (features above a certain
#' fraction of the Bernoulli variance) in 'bsseq' and returns their values.
#'
#' For DNA methylation, particularly when summarized across regions, we can do
#' better (a lot better) than MAD. Since we know:
#' max(SD(X_j)) if X_j ~ Beta(a, b) < max(SD(X_j)) if X_j ~ Bernoulli(a/(a+b))
#' for X with a known mean and standard deviation (SD), then we can solve for
#' (a+b) by MoM. We can then define the extremality by:
#' extremality = sd(X_j) / bernoulliSD(mean(X_j))
#'
#' @param bsseq  A bsseq object
#' @param r      Regions to consider - NULL covers all loci (DEFAULT: NULL)
#' @param k      How many rows/regions to return (DEFAULT: 500)
#'
#' @return       A GRanges object with methylation values sorted by extremality
#'
#' @import SummarizedExperiment
#'
#' @examples
#'
#'   shuf_bed <- system.file("extdata", "MCF7_Cunha_chr11p15_shuffled.bed.gz",
#'                           package="biscuiteer")
#'   orig_bed <- system.file("extdata", "MCF7_Cunha_chr11p15.bed.gz",
#'                           package="biscuiteer")
#'   shuf_vcf <- system.file("extdata",
#'                           "MCF7_Cunha_shuffled_header_only.vcf.gz",
#'                           package="biscuiteer")
#'   orig_vcf <- system.file("extdata",
#'                           "MCF7_Cunha_header_only.vcf.gz",
#'                           package="biscuiteer")
#'   bisc1 <- readBiscuit(BEDfile = shuf_bed, VCFfile = shuf_vcf,
#'                        merged = FALSE)
#'   bisc2 <- readBiscuit(BEDfile = orig_bed, VCFfile = orig_vcf,
#'                        merged = FALSE)
#'
#'   reg <- GRanges(seqnames = rep("chr11",5),
#'                  strand = rep("*",5),
#'                  ranges = IRanges(start = c(0,2.8e6,1.17e7,1.38e7,1.69e7),
#'                                   end= c(2.8e6,1.17e7,1.38e7,1.69e7,2.2e7))
#'                  )
#'
#'   comb <- unionize(bisc1, bisc2)
#'
#'   ext <- byExtremality(comb, r = reg)
#'
#' @export
#'
byExtremality <- function(bsseq,
                          r = NULL,
                          k = 500) {
 
  if (ncol(bsseq) == 1) {
    stop("No variance in a one sample object. Either use a data set ",
         "with multiple samples or combine bsseq objects with ",
         "unionize().")
  }

  if (is.null(r)) {
    message("No regions specified. This may melt your computer.")  
    byRegion <- FALSE
  } else { 
    byRegion <- TRUE 
  }

  if (nrow(bsseq) < k) {  
    message("Requested k (", k, ") exceeds the row count in bsseq ",
            "Adjusting k.") 
    k <- min(nrow(bsseq), k)
  }

  if (!is.null(r) & length(r) < k) {
    message("Requested k (", k, ") exceeds the region count. Adjusting k.") 
    k <- min(length(r), k)
  }

  if (byRegion) {
    bsseq <- sort(bsseq) 
    r <- sort(r) 
    m <- getMeth(bsseq, regions=r, type="raw", what="perRegion")
    rownames(m) <- as.character(r) 
  } else { 
    m <- getMeth(bsseq, type="raw")
    rownames(m) <- as.character(rowRanges(bsseq)) 
  }
  
  extr <- extremality(m)
  chosen <- rev(order(extr, na.last=FALSE))[seq_len(k)]
  res <- m[chosen, ]
  res <- cbind(res, extremality = extr[chosen])

  return(.makeGRangesFromMatrix(res))

}

# Matrix must have rownames from as.character(granges)
.makeGRangesFromMatrix <- function(matrix) {
  if (is.null(rownames(matrix))) {
    stop("Input matrix must have rownames corresponding to",
         "as.character(granges), where granges is a GRanges object")
  }
  # Start with the names of the GRanges from as.character(granges)
  df <- data.frame(grNames = rownames(matrix))
  # Split off the chromosome portion from ranges portion
  df <- data.frame(df, data.frame(do.call('rbind',
                                          strsplit(as.character(df$grNames),':',
                                                   fixed=TRUE))))
  # Make names easier to work with
  names(df) <- c("grNames","chr","range")
  # Split up the start and end parts from the ranges portion
  df <- data.frame(df, data.frame(do.call('rbind',
                                          strsplit(as.character(df$range),'-',
                                                   fixed=TRUE))))
  # Rename data.frame to straighforward names
  names(df) <- c("grNames","chr","range","start","end")
  df <- df[, c("chr","start","end")]

  df <- data.frame(df, matrix)

  return(makeGRangesFromDataFrame(df,keep.extra.columns=TRUE))
}
