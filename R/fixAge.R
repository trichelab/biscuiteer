#' Turn 'epigenetic clock' into actual age
#'
#' Uses Horvath-type 'epigenetic clock' raw output to project into actual ages
#'
#' The 'Epigenetic Clock' (Horvath 2012) and similar schemes use a number of 
#' CpG loci (or regions, or perhaps CpH loci -- it doesn't really matter what)
#' to estimate the chronological/biological age of samples from DNA methylation
#' with pre-trained feature weights (coefficients) for each region/locus. 
#' 
#' All of these type of clocks use a nonlinear output transformation which 
#' switches from an exponential growth model for children into a linear model
#' for adults, where `adult` is an arbitrary number (by default and custom,
#' that number is 21; elsewhere it can sometimes be seen as 20, but all known 
#' epi-age transformation functions quietly add 1 to the constant internally). 
#'
#' This function implements the above standard output transformation step.
#' 
#' @param x      Untransformed or raw prediction(s)
#' @param adult  Age of adulthood (DEFAULT: 21)
#' 
#' @return       Transformed prediction(s)
#'
#' @examples
#'
#'   clock <- getClock(genome="hg38")
#'   score <- clock$gr$score
#'
#'   fixAge(score)
#' 
#' @export
#'
fixAge <- function(x,
                   adult = 21) {
  return(ifelse(x < 0, adult*exp(x) - 1, adult*x + adult))
}
