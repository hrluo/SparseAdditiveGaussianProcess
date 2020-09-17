############################
# SAGP Model  - Kernel part#
############################

message('kernel.R module loaded!\n')

# Correlation functions  
#-----------------------
# From Pratols's STAT8810, 2017
# For now, we only use the gaussian kernel (powexp with alpha=2)

# powexp
# Power exponential, alpha=2 (default) is gaussian
powexp <-function(geoD,rho,alpha=2)
{
  rho=as.vector(rho)
  if(!is.vector(rho)) stop("non-vector rho!")
  if(any(rho<0)) stop("rho<0!")
  if(any(rho>1)) stop("rho>1!")
  if(any(alpha<1) || any(alpha>2)) stop("alpha out of bounds!")
  if(!is.list(geoD)) stop("wrong format for distance matrices list!")
  if(length(geoD)!=length(rho)) stop("rho vector doesn't match distance list")
  
  R=matrix(1,nrow=nrow(geoD$l1$m1),ncol=ncol(geoD$l1$m1))
  for(i in 1:length(rho))
    R=R*rho[i]^(geoD[[i]][[1]]^alpha)
  
  return(list(R=R))
}

# matern32
# Matern 3 2
matern32<-function(geoD,theta)
{
  theta=as.vector(theta)
  if(!is.vector(theta)) stop("non-vector theta!")
  if(any(theta<0)) stop("theta<0!")
  if(!is.list(geoD)) stop("wrong format for distance matrices list!")
  if(length(geoD)!=length(theta)) stop("theta vector doesn't match distance list")
  
  R=matrix(1,nrow=nrow(geoD$l1$m1),ncol=ncol(geoD$l1$m1))
  for(i in 1:length(theta))
  {
    D=(1+sqrt(3)*geoD[[i]][[1]]/theta[i])*exp(-sqrt(3)*geoD[[i]][[1]]/theta[i])
    R=R*D
  }
  
  return(list(R=R))
}

# matern52
# Matern 5 2
matern52<-function(geoD,theta)
{
  theta=as.vector(theta)
  if(!is.vector(theta)) stop("non-vector theta!")
  if(any(theta<0)) stop("theta<0!")
  if(!is.list(geoD)) stop("wrong format for distance matrices list!")
  if(length(geoD)!=length(theta)) stop("theta vector doesn't match distance list")
  
  R=matrix(1,nrow=nrow(geoD$l1$m1),ncol=ncol(geoD$l1$m1))
  for(i in 1:length(theta))
  {
    D=(1+sqrt(3)*geoD[[i]][[1]]/theta[i]+5*geoD[[i]][[1]]^2/(5*theta[i]^2))*exp(-sqrt(5)*geoD[[i]][[1]]/theta[i])
    R=R*D
  }
  
  return(list(R=R))
}

# wendland1
# Wendland 1
wendland1<-function(geoD,theta)
{
  theta=as.vector(theta)
  if(!is.vector(theta)) stop("non-vector theta!")
  if(any(theta<0)) stop("theta<0!")
  if(!is.list(geoD)) stop("wrong format for distance matrices list!")
  if(length(geoD)!=length(theta)) stop("theta vector doesn't match distance list")
  
  R=matrix(1,nrow=nrow(geoD$l1$m1),ncol=ncol(geoD$l1$m1))
  for(i in 1:length(theta))
  {
    D=(1-geoD[[i]][[1]]/theta[i])
    D[D<0]=0
    D=D^4*(1+4*geoD[[i]][[1]]/theta[i])
    R=R*D
  }
  
  return(list(R=R))	
}

# wendland2
# Wendland 2
wendland2<-function(geoD,theta)
{
  theta=as.vector(theta)
  if(!is.vector(theta)) stop("non-vector theta!")
  if(any(theta<0)) stop("theta<0!")
  if(!is.list(geoD)) stop("wrong format for distance matrices list!")
  if(length(geoD)!=length(theta)) stop("theta vector doesn't match distance list")
  
  R=matrix(1,nrow=nrow(geoD$l1$m1),ncol=ncol(geoD$l1$m1))
  for(i in 1:length(theta))
  {
    D=(1-geoD[[i]][[1]]/theta[i])
    D[D<0]=0
    D=D^6*(1+6*geoD[[i]][[1]]/theta[i]+35*geoD[[i]][[1]]^2/(3*theta[i]^2))
    R=R*D
  }
  
  return(list(R=R))	
}

# generalized.wendland
# Generalized Wendland 
generalized.wendland<-function(geoD,theta,kap)
{
  d=length(geoD)
  mu=(d+1)/2  # strictly speaking we need mu>=(d+1)/2 but we can change kappa also so
  # this is fair to assume.
  
  if(length(theta)>1) stop("theta is incorrect dimension\n")
  if(length(kap)>1) stop("kappa is incorrect dimensions\n")
  if(mu<((d+1)/2) ) stop("mu does not satisfy constraints\n")
  if(kap<0) stop("kappa > 0 required\n")
  
  D=matrix(0,nrow=nrow(geoD$l1$m1),ncol=ncol(geoD$l1$m1))
  for(i in 1:length(geoD))
    D=D+(geoD[[i]][[1]]^2)
  D=sqrt(D)
  
  if(kap==0) {
    # kap=0 is essential the Askey correlation
    D=D/theta
    R=1-D
    R[R<0]=0
    R=R^(mu+kap)
  }
  else
  {
    # library(fields) implements the general case
    R=fields::Wendland(D,theta=theta,dimension=d,k=kap)
  }
  
  rm(D)
  return(list(R=R))
}

# makeKernel
# Generate covariance matrix between X1 and X2
# kinfo is a list object that includes the parameters of the kernel (rho and eta for gaussian).
makeKernel<-function(X1,X2,kinfo){
  
  D <- ncol(X1)
  
  #Parameters for covariance matrix
  eta <- kinfo$eta
  rho <- kinfo$rho
  
  #geoD is list used in computing correlation matrices.
  #each element of the list is a matrix of differences and corresponds to one dimension.
  geoD <- list()
  for(d in 1:D) {
    geoD$ltemp <- list(mtemp=outer(as.vector(X1[,d]), as.vector(X2[,d]), "-"))
    names(geoD) <- paste0("l",1:d)
    names(geoD[[d]]) <- paste0("m",d)
  }
  
  #R is the correlation matrix. For now, only gaussian.
  R <- powexp(geoD, rho = rho, alpha = 2)
  
  #Covariance matrix =  1/precision * correlation mat
  K <- 1/eta * R$R
  
  return(K)
}
