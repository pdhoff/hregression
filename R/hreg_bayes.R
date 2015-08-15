#' Bayes Linear Regression Posterior Simulation
#'
#' This functions returns a simulation from the posterior distribution
#' of the vector of multivarate linear regression coefficients.
#'
#' @param y An \code{n} by \code{1} matrix.
#' @param X An \code{n} by \code{p} matrix.
#' @param b0 A \code{p} by \code{1} matrix.
#' @param V0 A \code{p} by \code{p} matrix.
#' @param ve An scalar.
#' @param Xty The value of \code{crossprod(X,y)}.
#' @param tXX The value of \code{crossprod(X)}.
#' @param iV0b0 The value of \code{ solve(V0)\%*\%b0 }.
#' @param iV0 The value of \code{ solve(V0) }.
#'
#' @return A \code{p} by \code{1} vector simulated 
#' from the posterior distribution.
#'
#' @author Peter Hoff
#'
#' @export
#'
#' @examples
#' p<-3 ; n<-10
#' b<-rnorm(p)
#' X<-rsan(c(n,p))
#' e<-rnorm(n)
#' y<-X%*%b+e
#'
#' rb_fc(y,X,rep(0,p),diag(nrow=p),1)
#'
#' mvreg_ols(y,X)
#'
rb_fc<-function(y,X,b0,V0,ve=1,
  tXy=crossprod(X,y),
  tXX=crossprod(X),
  iV0b0=solve(V0)%*%b0,
  iV0=solve(V0))
{

  Vb<-solve( iV0 + tXX/ve )
  eb<-Vb%*%( iV0b0 + tXy/ve )
  sb<-eb + t(chol(Vb))%*%rnorm(length(b0))
  # could possibly speed up with chol, chol2inv

  sb
}





#' Bayes Hierarchical Linear Regression Posterior Simulation
#'
#' This functions returns a simulation from the posterior distribution
#' of the matrix of linear regression coefficients.
#'
#' @param Y An \code{r} by \code{n} matrix
#' @param X An \code{r} by \code{p} by \code{n} array.
#' @param b0 A \code{p} vector.
#' @param V0 An \code{p} by \code{p} matrix.
#' @param ve A scalar. 
#'
#' @return An \code{r} by \code{p} matrix, simulated
#' from the posterior distribution.
#'
#' @author Peter Hoff
#' 
#' @export
#'
rB_fc<-function(Y,X,b0,V0,ve=1)
{
  iV0<-solve(V0) ; iV0b0<-iV0%*%b0

  B<-matrix(nrow=nrow(Y),ncol=dim(X)[2]) 
  for(i in 1:dim(Y)[1])
  {
    B[i,]<-rb_fc(c(Y[i,]),t(array(X[i,,],dim=dim(X)[-1])), 
                 b0,ve=ve,iV0b0=iV0b0,iV0=iV0)
  }

  B
}



#' Array matrix slice product
#' 
#' Multiply matrices and vectors that 
#' correspond to the same slices as 
#' an array and matrix, respectively. 
#'
#' @param X An \code{r} by \code{p} by \code{n} array.
#' @param B An \code{r} by \code{p} matrix.
#' 
#' @return An \code{r} by \code{n} matrix, obtained by
#' matrix multiplying each row of \code{B} by the 
#' matrix obtained from the \code{r}th first-mode
#' slice of \code{X}. 
#' 
#' @author Peter Hoff
#' 
#' @export
amsprod<-function(X,B)
{
  apply( sweep(X,c(1,2),B,"*") ,c(1,3),sum)
}



#' Bayes Hierarchical Linear Regression Posterior Simulation
#'
#' This functions returns a simulation from the posterior distribution
#' of the shared error variance.
#'
#' @param Y An \code{r} by \code{n} matrix.
#' @param X An \code{r} by \code{p} by \code{n} array.
#' @param B An \code{r} by \code{p} matrix.
#'
#' @return A scalar simulated 
#' from the posterior distribution.
#'
#' @author Peter Hoff
#'
#' @export
#'
rve_fc<-function(Y,X,B,ss0=1,nu0=1)
{ 
  XB<-amsprod(X,B) 

  ssd<-sum((Y-XB)^2)

  nud<-prod(dim(Y))  

  1/rgamma(1,(nu0+nud)/2,(ss0 + ssd)/2 ) 
}


#' Bayes Hierarchical Linear Regression Posterior Simulation
#'
#' This functions returns a simulation from the posterior distribution
#' of the across-group covariance of the regression vectors. 
#'
#' @param B An \code{r} by \code{p} array.
#' @param b0 A \code{p} vector. 
#'
#' @return An \code{p} by \code{p} matrix simulated
#' from the posterior distribution.
#'
#' @author Peter Hoff
#'
#' @export
#'
rV0_fc<-function(B,b0,S0=diag(nrow=length(b0)),nu0=length(b0)+1)
{ 
  Sd<-crossprod( sweep(B,2,b0,"-") ) 

  nud<-nrow(B) 

  solve( rwish( solve(S0+Sd), nu0+nud ) )
}


#' Bayes Hierarchical Linear Regression Posterior Simulation
#'
#' This functions returns a simulation from the posterior distribution
#' of the across-group mean of the regression vectors. 
#'
#' @param B An \code{r} by \code{p} matrix
#' @param V0 An \code{mp} by \code{mp} matrix. 
#' @param iV0 The inverse of \code{V0}. 
#'
#' @return A \code{p} vector simulated
#' from the posterior distribution.
#'
#' @author Peter Hoff
#'
#' @export
#'
rb0_fc<-function(B,V0,iV0=solve(V0),
                 eb0=rep(0,ncol(B)),iVb0=diag(nrow=ncol(B))/nrow(B) )  
{ 

  Vb0p<-solve( iVb0 + iV0*nrow(B) )
  eb0p<-Vb0p%*%( iVb0%*%eb0 + iV0%*%apply(B,2,sum) ) 

  eb0p+t(chol(Vb0p))%*%rnorm(ncol(B)) 
}


