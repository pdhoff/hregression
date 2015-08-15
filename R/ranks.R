#' Kendall's tau measure of association 
#' 
#' This function provides a Monte Carlo approximation to 
#' Kendall's tau measure of association. 
#' 
#' @param x a vector
#' @param y a vector
#' @param nmc an integer number of Monte Carlo simulations
#' 
#' @return An array of dimension \code{dim}. 
#'
#' @author Peter Hoff
#' 
#' @export 
#'  
#' @examples
#' rsan(c(5,4,3))
#' 
kendalltau<-function(x,y,nmc=1e5)
{
  i<-sample(length(y),nmc,replace=TRUE)
  j<-sample(length(y),nmc,replace=TRUE)
  s1<-cbind(x[i],y[i])
  s2<-cbind(x[j],y[j])
  cp<-sum( (s1[,1]-s2[,1])*(s1[,2]-s2[,2])  >0 ,na.rm=TRUE)
  dp<-sum( (s1[,1]-s2[,1])*(s1[,2]-s2[,2])  <0 ,na.rm=TRUE)

  (cp-dp)/sum(cp+dp)
}


#' Normal Scores
#' 
#' This function applies a quantile-quantile transformation to 
#' the data, resulting in a distribution that is approximately normal
#' but has the same ranks as the original data
#' 
#' @param y A vector. 
#' @param ties.method The option \code{ties.method} in the \code{rank} function.
#' 
#' @return A vector of the same length as \code{y}. 
#' 
#' @author Peter Hoff
#' 
#' @export
#'
#' @examples
#' y<-rexp(100)
#' z<-zscores(y) 
#' par(mfrow=c(1,3))
#' hist(y) 
#' hist(z)
#' plot(y,z) 
zscores<-function(y,ties.method="average")
{
  z<-qnorm(rank(y,na.last="keep",ties.method=ties.method)/(sum(!is.na(y))+1) )

  names(z)<-names(y) 

  m<-dim(y)
  if(length(m)==2){ z<-matrix(z,nrow=m[1],ncol=m[2]) ; dimnames(z)<-dimnames(y) }
  if(length(m)>=3){ z<-array(z,dim=m) ; dimnames(z)<-dimnames(y) } 

  z 
}





#' Gibbs sampling for the Gaussian rank likelihood
#'
#' This function simulates a latent Gaussian object
#' \code{z} conditional on a partial ordering induced by 
#' an observed data object \code{y}. 
#'
#' @param y A vector, matrix or array. 
#' @param z A vector, matrix or array of the same dimension as \code{y}. 
#' @param ez A vector, matrix or array of the same dimension as \code{y}. 
#' @param sz A scalar, the common standard deviation of \code{z}. 
#'
#' @return A vector, matrix or array of the same dimension as \code{y}. 
#'
#' @author Peter Hoff
#'
#' @export
#'
rz_fc<-function(y,z,ez,sz=1)
{ 
  m<-dim(y) ; y<-c(y) ; z<-c(z) ; ez<-c(ez) 

  r<-match(y,sort(unique(y))) 
  ikm<-numeric(0)
  ik0<-which(r==1)
  for(k in 1:max(r,na.rm=TRUE))
  {
    ikp<-which(r==k+1)

    lb<-suppressWarnings( max( z[ikm], na.rm=TRUE ) )
    ub<-suppressWarnings( min( z[ikp], na.rm=TRUE ) )

    mk<-ez[ik0] 
   
    zk<-mk+sz*qnorm(runif(length(ik0), pnorm((lb-mk)/sz), pnorm((ub-mk)/sz) ))

    # check for numerical error
    isinf<-which(abs(zk)==Inf | zk<lb | zk>ub)
    if(length(isinf)>0)
    {
      zc<-z[ik0][isinf]
      eps<-min(ub-lb,1)
      zp<-rnorm(length(zc),zc,eps/4)

      lu<-log(runif(length(zc)))
      lr<-dnorm(zp,mk[isinf],sz,log=TRUE) - 
          dnorm(zc,mk[isinf],sz,log=TRUE) +
          log( lb<zp & zp<ub)

      zc[lu<lr]<-zp[lu<lr]
      zk[isinf]<-zc
    } 

    z[ik0]<-zk
    ikm<-ik0 ; ik0<-ikp
  }


  # unconstrained if y is missing
  z[is.na(r)]<-rnorm(sum(is.na(r)),ez[is.na(r)],sz) 


  # global scale update
  sc<-sqrt(sum(z^2))
  u<-z/sc ; a<-length(z)-1 ; b<-sum(u*ez)/sz^2 ; c<-1/sz^2

  es<-(b+sqrt(b^2+4*a*c))/(2*c) 
  ss<-sqrt( 1/(a/es^2 + c) )      
  sp<-rnorm(1,es,ss)  

  lhr<- a*log(sp/sc) + b*(sp-sc) - c*.5*(sp^2-sc^2) +
        dnorm(sc,es,ss,log=TRUE) - dnorm(sp,es,ss,log=TRUE ) + log(sp>0)
  if(log(runif(1))<lhr) { sc<-sp } 
  
  z<-u*sc 


  # global location update
  z<-(z-mean(z)) + rnorm(1, mean(ez), 1/length(ez)) 


  # reshape
  if(length(m)==2){ z<-matrix(z,nrow=m[1],ncol=m[2]) } 
  if(length(m)>=3){ z<-array(z,dim=m) }   
  z

}



