#' a simple local wrapper around HDF5Array::saveHDF5SummarizedExperiment
#' 
#' @param   x     the SummarizedExperiment-derived object to save as HDF5 
#' @param   dir   path to store the HDF5 file "assays.h5" (default is "HDF5") 
#' @param   ...   arguments to pass HDF5Array::saveHDF5SummarizedExperiment 
#' 
#' @return        invisibly, the HDF5-backed version of the object being saved
#' 
#' @import SummarizedExperiment
#' @import HDF5Array
#' 
#' @export
saveHDF5 <- function(x, dir="HDF5", ...) {
  # lop off the filename if it is supplied
  dir <- base::sub("assays.h5$", "", dir)
  saveHDF5SummarizedExperiment(x, dir=dir, ...) 
}
