#' helper function to play nice with BSseq's expectations about the M matrix 
#' 
#' @param x       the matrix-like object containing NAs (holes) to fix
#' @param y       the value with which to fill the holes (default is 0)
#' @param sparse  make the result Matrix-backed? (NULL; if beneficial)
#'
#' @return    x, but with holes filled, and as a matrix
#'
#' @export
fixNAs <- function(x, y=0, sparse=NULL) { 
  if (is.null(sparse) && median(x, na.rm=TRUE) < 0.01) {
    sparse <- TRUE
  }
  if (sparse) {
    x <- sparseMatrix(x)
  } else {
    x <- as.matrix(x)
  }
  x[is.na(x)] <- y
  return(x)
}
