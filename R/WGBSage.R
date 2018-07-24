#' guess a specimen's age from a modified Horvath algorithm
#' 
#' @param   x     a BSseq object
#'  
#' @return  age   age estimate for each column of x 
#'
#' @export
WGBSage <- function(x) { 
 
  stopifnot(unique(genome(x))[1] %in% c("hg19","GRCh37"))
  data("Horvath_EPIC.hg19", package="biscuitEater") 
  
  forWGBSage <- getMeth(x, regions=Horvath_EPIC.hg19, 
                        type="raw", what="perRegion")
  message("Still not quite done") 

  return(forWGBSage)
    
}
