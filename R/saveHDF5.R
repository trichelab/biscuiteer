#' Save HDF5-backed SummarizedExperiment
#'
#' Local wrapper around HDF5Array::saveHDF5SummarizedExperiment: NOTE: CURRENTLY
#' BISCUITEER DOES NOT HAVE HDF5-BACKING FUNCTIONALITY. USE AT YOUR OWN RISE!!!!
#'
#' @param x    The SummarizedExperiment-derived object to save as HDF5
#' @param dir  Path to store the HDF5 file "assays.h5" (DEFAULT: "HDF5")
#' @param ...  Arguments to pass to HDF5Array::saveHDF5SummarizedExperiment
#'
#' @return     Invisible. The HDF5-backed version of the object being saved
#'
#' @importFrom HDF5Array saveHDF5SummarizedExperiment
#'
#' @examples
#'
#' @export
#'
saveHDF5 <- function(x,
                     dir = "HDF5",
                     ...) {
  # lop off the filename if it is supplied
  dir <- base::sub("assays.h5$", "", dir)
  HDF5Array::saveHDF5SummarizedExperiment(x, dir=dir, ...) 
}
