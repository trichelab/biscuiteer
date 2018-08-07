#' guess a specimen's age from a modified Horvath algorithm
#' 
#' @param   x     a BSseq object
#'  
#' @return  age   age estimate for each column of x 
#'
#' @export
WGBSage <- function(x) { 
 
  stopifnot(unique(genome(x))[1] %in% c("hg19","GRCh37"))
  data("Horvath_EPIC.hg19", package="biscuiteer") 
  forWGBSage <- getMeth(x, regions=Horvath_EPIC, type="raw", what="perRegion")
  rownames(forWGBSage) <- Horvath_EPIC$name
  message("Calibration of weights still needs tuning. Returning raw 5mC%...") 
  return(forWGBSage)
    
}
