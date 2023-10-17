### deshrinking ridge regression for multi-locus analyses
#' @param x is a 1 vector, which is a covariate
#' @param y is a matrix of traits
#' @param z is a matrix of genotypes
#' @param qq is the eigenvalue decomposition of matrix z 

ridge<-function(x,y,z,qq){
  
  u<-qq$vectors
  d<-qq$values
  q<-ncol(x)
  xu<-t(u)%*%x
  yu<-t(u)%*%y
  zu<-t(u)%*%z
  n<-nrow(z)
  m<-ncol(z)
  
  loglike<-function(theta){
    lambda<-exp(theta)
    h<-sqrt(1/(lambda*d+1))
    x1<-xu*h
    y1<-yu*h
    xHx<-t(x1)%*%x1
    xHy<-t(x1)%*%y1
    b<-solve(xHx,xHy)
    ru<-yu-xu%*%b
    r1<-ru*h
    rHr<-t(r1)%*%r1
    v<-drop(rHr)/(n-q)
    d1<-sum(log(lambda*d+1))
    d2<-log(det(xHx))
    d3<-(n-q)*log(v)
    loglike<- 0.5*(d1+d2+d3)
    return(loglike)
  }
  
  fixed<-function(theta){
    lambda<-exp(theta)
    h<-sqrt(1/(lambda*d+1))
    x1<-xu*h
    y1<-yu*h
    z1<-zu*h
    xHx<-t(x1)%*%x1
    xHy<-t(x1)%*%y1
    zHz<-apply(z1^2,2,sum)
    
    c11<-solve(xHx)
    b<-c11%*%xHy
    ru<-yu-xu%*%b
    r1<-ru*h
    rHr<-t(r1)%*%r1
    ve<-drop(rHr)/(n-q)
    vg<-lambda*ve
    varb<-diag(c11*ve)
    gamma<-lambda*t(z1)%*%r1
    d22<-vg*zHz*lambda
    
    c22<-rep(vg,m)-d22
    
    vgk1<-c22
    
    vgk2<-d22
    
    df<-zHz*lambda
    
    covparm<-data.frame(vg,ve,lambda)
    estimate<-data.frame(b,varb)
    blup<-data.frame(gamma,vgk1,vgk2,df)
    
    return(list(covparm,estimate,blup))
  }
  
  theta0<-1
  parm<-optim(par=theta0,fn=loglike,method="L-BFGS-B",lower=-Inf,upper=Inf)
  bb<-fixed(theta=parm$par)
  return(bb)
}
