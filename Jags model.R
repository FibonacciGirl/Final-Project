install.packages("runjags")

library(runjags)


data<-read.csv('TD')
category<-data$category
data$category<-NULL
n.categories<-length(unique(category))
Y<-data.matrix(data)



data.model<-list(Y=Y, Sigma0<-cov(Y), n = dim(Y)[1], n.dimensions = dim(Y)[2], n.categories=n.categories)
result<- run.jags('JAGS.txt', data=data.model, monitor=c('theta','W','mu','sigma'), n.chains=3, burnin=1000, sample=10000)
