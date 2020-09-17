source("Core/sagpfuns.r")
source("Core/kernel.R")

########
# Data #
########

#Generate data to fit SAGP model. Uncomment the portion of code 
#to generate the desired type of data.
#Parameters, determine the location where new predictions lie.
edgeGridX <- seq(0,1,length = 100)
edgeGridY <- seq(0,1,length = 100)
# In terms of UTM coordinates.
Xnew <- as.matrix(expand.grid(edgeGridX,edgeGridY))
truthFUN<-function(X1,X2){X1*sin(X2)}
X<-expand.grid(seq(-10,10,1),seq(-10,10,1))
X<-as.matrix(X)
truth<-c()
Y<-c()
for(j in 1:nrow(X)){
  truth[j]=truthFUN(X[j,1],X[j,2])
  Y[j]=truth[j]+rnorm(1)
}
prep<-read.table('test.txt')
Y<-prep[,1]
Y<-as.vector(Y)
X<-prep[,2:3]
X<-as.matrix(X)


#range01 is a force-scaled function that stipulates the scale to be on [0,1].
X1t_dat<-min(X[,1])
X2t_dat<-min(X[,2])
X1m_dat<-max(X[,1])-min(X[,1])
X2m_dat<-max(X[,2])-min(X[,2])
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
X[,1]<-range01(X[,1])
X[,2]<-range01(X[,2])

X1new <- edgeGridX
X2new <- edgeGridY
#Make marginal mean and variance 0 and 1
meanY <- mean(Y)
varY <- as.numeric(var(Y))
Y = (Y - meanY)/sqrt(varY)
para=list(Yt=meanY,Ym=sqrt(varY),
          X1t=X1t_dat,X1m=X1m_dat,
          X2t=X2t_dat,X2m=X2m_dat)



##########
# Priors #
##########

#Hyperparameters for priors.

pinfo=list(m=5, #number of pseudo-inputs
           a=100, b=1, # hyp. for lambda (precision of noise) - gamma prior
           r=.1, #Parameter that controls prior of eta
           logRhoFirstLayer=-1, # logarithm (base 10) of rho in first layer
           logRhoLastLayer=-50  # logarithm (base 10) of rho in last layer
)
########
# MCMC #
########
#Parameters of the MCMC sampler.
minfo=list(Nmcmc=10, #num saved samples
           burnIn=10, #num discarded samples for burn in
           semiWidthRho=0.05, #semiwidth of the proposal uniform distribution in MH step for rho (tuned adaptively)
           semiWidthEta=0.05, #semiwidth of the proposal uniform distribution in MH step for eta (tuned adaptively)
           seed=1234) #seed of the sampler, for reproducibility

tinfo <- generateTree(list(c(0,1),c(0,1)),L = 2)

#############
# Fit Model #
#############
starting <- Sys.time()
fit <- sagp(Y, X, pinfo, minfo, tinfo, Xnew=X, sampleEta=T, sampleRho=F, storePredX=T)
ending <- Sys.time()

############
# Save PDF #
############
library(ggplot2)
ggplot()+geom_point(aes(x=X[,1],y=X[,2],col=Y),size=5)+scale_colour_viridis_c(limits = c(0, 1))
ggplot()+geom_point(aes(x=fit$Xnew[,1],y=fit$Xnew[,2],col=rowMeans(fit$mu.new)),size=5)+scale_colour_viridis_c(limits = c(-1, 1))
ggplot()+geom_point(aes(x=fit$Xnew[,1],y=fit$Xnew[,2],col=rowMeans(fit$Yhat)),size=5)+scale_colour_viridis_c(limits = c(-1, 1))





