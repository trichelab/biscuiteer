#' helper function to play nice with BSseq's expectations about the M matrix 
#' 
#' @param x       the matrix-like object containing NAs (holes) to fix
#' @param y       the value with which to fill the holes (default is 0)
#' @param sparse  make the result Matrix-backed? (TRUE)
#'
#' @return x, but with holes filled, and as a (possibly sparse) matrix
#'
#' @importFrom Matrix Matrix
#' 
#' @export
fixNAs <- function(x, y=0, sparse=TRUE) { 
  for (i in seq_len(ncol(x))) x[which(is.na(x[,i])), i] <- y
  x <- as.matrix(x)
  if (sparse) x <- Matrix(x)
  return(x)
}
