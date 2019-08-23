#' a simple local wrapper around HDF5Array::loadHDF5SummarizedExperiment
#' 
#' @param   dir   path to the file assays.h5 (and se.rds); default is "HDF5"
#' 
#' @return        a SummarizedExperiment-derived object (such as a BSseq object)
#' 
#' @import HDF5Array
#' 
#' @export
loadHDF5 <- function(dir="HDF5") {
  # lop off the filename if it is supplied
  dir <- base::sub("assays.h5$", "", dir)
  loadHDF5SummarizedExperiment(dir=dir)
}
