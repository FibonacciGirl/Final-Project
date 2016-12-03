#install.packages("mvtnorm")
#install.packages("MCMCpack")
library(mvtnorm)
library(MCMCpack)
library(coda)


##credit Ravi Varadhan rvaradhan at jhmi.edu ####
Posdef <- function (r, ev = runif(r, 0, 1)) 
{
  U <- matrix(ncol=r, rnorm(r^2))
  decomp <- qr(U)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  H <- t(O) %*% diag(ev) %*% O
  return(H)
}  

Posdef(2)
######################

n.categories<-2 # run the loop this many times to get prob for each category
n.dimensions<-2 #number of dimensions in which the stimulus is observed
training.data<-
test.data<-



Posdef(n.dimensions)
  
  
  #latent class stuff#



category.probability<-matrix(ncol=n.categries,nrow=n.run)

MCMC<-function(Y,ybar,mu0,Lambda0,nu0,sigma0,alpha){
  n.categories<-2
  n<-dim(Y)[2]
  n.dimensions<-dim(Y)[1]
  
  
  ##MCMC setup##
  nrun=10000
  thin=
    burn=


  ##begin sampler##
  for(iter in 2:nrun){
    print(iter)
    ##update Z##
    Z.jump<-array(dim=c(n,n.categories))
    
        for(s in 1:n){
          X<-numeric(n.categories)
          for(j in 1:n.categories){
            X[j]<-rgamma(1,alpha[j])
          }
          X.sum<-sum(X)
          P<-numeric(n.categories)
          for(j in 1:n.categories){
            P[j]<-X[j]/X.sum
          }
          Z.jump[s,]<-rmultinom(1,1,prob=P)
        }
    loglik.current<-numeric(n)
    loglik.jump<-numeric(n)
    prior.current<-numeric(n)
    prior.jump<-numeric(n)
    
    for(s in 1:n){
      j.jump<-which(Z.jump[s,]==1)
      j.current<-which(Z[(iter-1),s,]==1)
      loglik.current[s]<-log(dmvnorm(Y[,s],Mu[(iter-1),j.current,],Sigma[(iter-1),j.current,,]))
      loglik.jump[s]<-log(dmvnorm(Y[,s],Mu[(iter-1),j.jump,],Sigma[(iter-1),j.jump,,]))
      prior.current[s]<-log(dmultinom(Z[(iter-1),s,],1,theta[(iter-1),]))
      prior.jump[s]<-log(dmultinom(Z.jump[s,],1,theta[(iter-1),]))
    }
    
    prob.jump<-sum(loglik.jump)+sum(prior.jump)
    prob.current<-sum(loglik.current)+sum(prior.current)
    
    accept.prob<- min(1,exp(prob.jump)/exp(prob.current))
    print(accept.prob)
    
    accept<-rbinom(1,1,accept.prob)
    if(accept == 1){
      Z[iter,,]<-Z.jump
    }
    if(accept == 0){
      Z[iter,,]<-Z[(iter-1),,]
    }

    ##update theta##
    Z.sum<-numeric(n.categories)
    
    for(j in 1:n.categories){
      Z.sum[j]<-sum(Z[iter,,j])
    }
    
    theta[iter,]<-theta[(iter-1),]+Z.sum
    
    print(theta[iter,])
    ##update mu##
    for(j in 1:n.categories){
      Lambdan<-solve(solve(Lambda0[j,,])+n*solve(Sigma[(iter-1),j,,]))
      mun<-Lambdan%*%((solve(Lambda0[j,,])%*%mu0[j,])+(n*solve(Sigma[(iter-1),j,,])%*%ybar[j,]))
      Mu[iter,j,]<-rmvnorm(1,mun,Lambdan)
    }
    

    
    ##update sigma##
    for(j in 1:n.categories){
      nun<-nu0 + n
      Sn<-sigma0[j,,]+(ybar[j,]-Mu[(iter-1),j,])%*%(t(ybar[j,]-Mu[(iter-1),j,]))

      sig.inv<-rWishart(1,nun,solve(Sn))
      Sigma[iter,j,,]<-solve(matrix(sig.inv,nrow=n.dimensions,ncol=n.dimensions))
    }
  }
  
  ##burn and thin##
  
  ##return##
}
##data##
data<-read.csv('TD')
category<-data$category
data$category<-NULL
n.categories<-length(unique(category))
NOPE<-t(data.matrix(data))
Y<-NOPE[,1:4]

n.categories<-2
n<-dim(Y)[2]
n.dimensions<-dim(Y)[1]

ybar<-array(dim= c(n.categories,n.dimensions))

for(j in 1:n.categories){
  ybar[j,]<-mean(Y[,which(category[1:4]==j)])
}

#choose fixed hyper paratments#
mu0<-array(dim=c(n.categories,n.dimensions))
mu0<-rmvnorm(n.categories,rep(0,n.dimensions),Posdef(n.dimensions))
Lambda0<-array(dim=c(n.categories,n.dimensions,n.dimensions))
for(j in 1:n.categories){
  Lambda0[j,,]<-rWishart(1,n.dimensions,Posdef(n.dimensions))
  
}

nu0<-numeric(n.categories)
nu0<-rep(n.dimensions,n.categories)
sigma0<-array(dim=c(n.categories,n.dimensions,n.dimensions))
for(j in 1:n.categories){
  sigma0[j,,]<-rWishart(1,n.dimensions,Posdef(n.dimensions))
}
alpha<-rep(1,n.categories)


##set initial values##
Z<-array(dim=c(nrun,n,n.categories))
theta<-array(dim=c(nrun,n.categories))
Mu<-array(dim=c(nrun,n.categories,n.dimensions))
Sigma<-array(dim=c(nrun,n.categories,n.dimensions,n.dimensions))

for(s in 1:n){
  Z[1,s,]<-rmultinom(1,1,rep(1/n.categories,n.categories))
}
theta[1,]<-rep(1/n.categories,n.categories)
for(j in 1:n.categories){
  Mu[1,j,]<-rmvnorm(1,mu0[j,],Lambda0[j,,])

  Sigma.inv<-rWishart(1,n.dimensions,sigma0[j,,])
  Sigma[1,j,,]<-solve(matrix(Sigma.inv,nrow=n.dimensions,ncol=n.dimensions))
}



MCMC(Y,ybar,mu0,Lambda0,nu0,sigma0,alpha)
  

