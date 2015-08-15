#' Wishart simulation
#'
#' Simulate a Wishart-distributed random matrix
#' using Bartletts decomposition, as described in Everson and
#' Morris (2000).
#'
#' @param S0 A positive definite matrix.
#' @param nu A positive scalar
#'
#' @return A square matrix. 
#' 
#' @examples
#' # simulate several matrices and compute the mean
#' SS<-matrix(0,5,5)
#' for(s in 1:1000) { SS<-SS+rwish(diag(5),3) }
#' SS/s
#'
#' @author Peter Hoff
#'
#' @export
rwish<-function(S0,nu=dim(as.matrix(S0))[1]+1)
{
  S0<-as.matrix(S0)
  S0h<-eigen(S0,symmetric=TRUE)
  S0h<-S0h$vec%*%diag(sqrt(S0h$val),nrow=length(S0h$val))%*%t(S0h$vec)

  p<-dim(S0)[1]
  T<-matrix(0,p,p)
  T[lower.tri(T)]<-rnorm( p*(p-1)/2)
  diag(T)<-sqrt( rgamma(p,(nu-(1:p)+1)/2,1/2) )
  S0h%*%T%*%t(T)%*%S0h
}

