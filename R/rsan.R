#' Random Standard Normal Array
#' 
#' This functions generates an array of dimension \code{dim}
#' filled with iid standard normal random variables. 
#' 
#' @param dim A vector of positive integers. 
#' 
#' @return An array of dimension \code{dim}. 
#' 
#' @examples
#' rsan(c(5,4,3))
#' 
#' @author Peter Hoff
#'
#' @export
rsan<-function(dim)
{
  array(rnorm(prod(dim)),dim)
}


