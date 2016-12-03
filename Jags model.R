install.packages("runjags")

library(runjags)
library(coda)

data<-read.csv('TD')
category<-data$category
data$category<-NULL
n.categories<-length(unique(category))
Y<-t(data.matrix(data))

dim(Y)
dim(L0)
mu0= matrix(c(.5,.5),nrow=2,ncol=1)
L0= matrix(c(0.4328403, 0.0406657, 0.0406657, 0.4222715),nrow=2,ncol=2)
nu0= 4


data.model<-list(Y=Y, S0=matrix(c(0.4328403, 0.0406657, 0.0406657, 0.4222715),nrow=2,ncol=2), nu0= 4, L0=matrix(c(0.4328403, 0.0406657, 0.0406657, 0.4222715),nrow=2,ncol=2),mu0= matrix(c(.5,.5),nrow=2,ncol=1), n = dim(Y)[1], p = dim(Y)[2], m=n.categories)
result<- run.jags('model.txt', data=data.model, monitor=c('theta','W','mu','sigma'), n.chains=3, burnin=1000, sample=10000)
