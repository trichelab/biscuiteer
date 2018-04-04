#' helper function to play nice with BSseq's expectations about the M matrix 
#' 
#' @param x       the matrix-like object containing NAs (holes) to fix
#' @param y       the value with which to fill the holes (default is 0)
#' @param sparse  make the result Matrix-backed? (FALSE)
#'
#' @return x, but with holes filled, and as a (possibly sparse) matrix
#'
#' @export
fixNAs <- function(x, y=0, sparse=FALSE) { 
  if (is.null(sparse) | sparse == FALSE) {
    x <- as.matrix(x)
  } else {
    x <- sparseMatrix(x)
  }
  x[is.na(x)] <- y
  return(x)
}
