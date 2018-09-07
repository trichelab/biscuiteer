#' choose loci or features by extremality (fraction of Bernoulli variance)
#'
#' For DNA methylation, particularly summarized across regions,
#' we can do better (a lot better) than MAD.  Since we know that 
#'
#' max(SD(X_j)) if X_j ~ Beta(a, b) < max(SD(X_j)) if X_j ~ Bernoulli(a/(a+b))
#'
#' for X having a known mean and SD, hence solvable for a + b by MoM, define
#'
#' extremality = sd(X_j) / bernoulliSD(mean(X_j))
#'
#' This function finds the k most extremal features in x & returns their values.
#' For count-based proportions, 
#'
#' @param     x           a bsseq object 
#' @param     r           regions to consider (default is NULL, all loci... !)
#' @param     k           how many rows/regions to return (500)
#' 
#' @return    a matrix-like object with methylation values, by extremality 
#' 
#' @export
byExtremality <- function(x, r=NULL, k=500) {
 
  if (is.null(r)) {
    message("No regions specified. This may melt your computer.")  
    byRegion <- FALSE
  } else { 
    byRegion <- TRUE 
  }

  if (nrow(x) < k) {  
    message("Requested k (", k, ") exceeds the row count in x. Adjusting k.") 
    k <- min(nrow(x), k)
  }

  if (!is.null(r) & length(r) < k) {
    message("Requested k (", k, ") exceeds the region count. Adjusting k.") 
    k <- min(length(r), k)
  }

  if (byRegion) {
    x <- sort(x) 
    r <- sort(r) 
    m <- getMeth(x, regions=r, type="raw", what="perRegion")
    rownames(m) <- as.character(r) 
  } else { 
    m <- getMeth(x, type="raw")
    rownames(m) <- as.character(rowRanges(x)) 
  }
  
  extr <- extremality(m)
  chosen <- rev(order(extr, na.last=FALSE))[seq_len(k)]
  res <- m[chosen, ]
  attr(res, "extremality") <- extr[chosen]
  return(res)

}
