#' Load HDF5-backed SummarizedExperiment
#'
#' Local wrapper around HDF5Array::loadHDF5SummarizedExperiment. NOTE: CURRENTLY
#' BISCUITEER DOES NOT HAVE HDF5-BACKING FUNCTIONALITY. USE AT YOUR OWN RISK!!!!  
#'
#' @param dir  Path to the assays.h5 (and se.rds) files (DEFAULT: "HDF5")
#'
#' @return     A SummarizedExperiment-derived object (list a bsseq object)
#'
#' @importFrom HDF5Array loadHDF5SummarizedExperiment
#'
#' @seealso HDF5Array::loadHDF5SummarizedExperiment
#'
#' @examples
#'
#' @export
#'
loadHDF5 <- function(dir = "HDF5") {
  # lop off the filename if it is supplied
  dir <- base::sub("assays.h5$", "", dir)
  HDF5Array::loadHDF5SummarizedExperiment(dir=dir)
}
