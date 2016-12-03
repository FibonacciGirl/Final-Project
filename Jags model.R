install.packages("rjags")

library(rjags)
library(runjags)
library(coda)

data<-read.csv('TD')
category<-data$category
data$category<-NULL
n.categories<-length(unique(category))
E<-t(data.matrix(data))
dim(E)

dim(L0)
mu0= matrix(c(.5,.5),nrow=2,ncol=1)
L0= matrix(c(0.4328403, 0.0406657, 0.0406657, 0.4222715),nrow=2,ncol=2)
W0= matrix(c(0.4328403, 0.0406657, 0.0406657, 0.4222715),nrow=2,ncol=2)


data.model<-list(E=E, S0=matrix(c(0.4328403, 0.0406657, 0.0406657, 0.4222715),nrow=2,ncol=2),W0=matrix(c(0.4328403, 0.0406657, 0.0406657, 0.4222715),nrow=2,ncol=2), L0=matrix(c(0.4328403, 0.0406657, 0.0406657, 0.4222715),nrow=2,ncol=2),mu0= matrix(c(.5,.5),nrow=2,ncol=1), n = dim(E)[2], p = dim(E)[1], m=n.categories)
result<- run.jags('model.txt', data=data.model, monitor=c('theta','W','mu','sigma'), n.chains=3, burnin=1000, sample=10000)
