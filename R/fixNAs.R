#' Replace NAs with another value
#'
#' Useful for coercing matrices into how bsseq is expecting the M matrix to be.
#'
#' @param x       The matrix-like object containing NAs to fix
#' @param y       The value to replace the NAs with (DEFAULT: 0)
#' @param sparse  Make the result a Matrix object? (DEFAULT: FALSE)
#'
#' @return        x with no NAs (possibly a sparse Matrix)
#'
#' @importFrom Matrix Matrix
#' 
#' @examples
#'
#'   nom <- c(rep(c(1,4,NA,9,NA,NA,7,NA), 5))
#'   fixNAs(nom)
#'
#' @export
#'
fixNAs <- function(x,
                   y = 0,
                   sparse = FALSE) { 
  x <- as(x, ifelse(sparse, "Matrix", "matrix"))
  x[is.na(x)] <- y
  return(x)
}
