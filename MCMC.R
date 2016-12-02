##install.packages("mvtnorm")
##install.packages("MCMCpack")
library(mvtnorm)
library(MCMCpack)
library(coda)


##credit Ravi Varadhan rvaradhan at jhmi.edu ####
Posdef <- function (n, ev = runif(n, 0, 1)) 
{
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  H <- t(O) %*% diag(ev) %*% O
  return(H)
}  
######################

n.categories<-4 # run the loop this many times to get prob for each category
n.dimensions<-3 #number of dimensions in which the stimulus is observed
training.data<-
test.data<-



Posdef(n.dimensions)
  
  
  #latent class stuff#

##data##
data<-read.csv('TD')
category<-data$category
data$category<-NULL
n.categories<-length(unique(category))
Y<-data.matrix(data)


##set initial values##
mu0<-array(dim=c(n.categories,n.dimensions))
mu0<-rmvnorm(n.categories,rep(0,n.dimensions),Posdef(n.dimensions))
Lambda0<-array(dim=c(n.categories,n.dimensions,n.dimensions))
for(j in 1:n.categories){
  Lambda0[j,,]<-rWishart(1,n.dimensions*n.dimensions,Posdef(n.dimensions))
  
}
nu0<-numeric(n.categories)
nu0<-rgamma(n.categories,1)
sigma0<-array(dim=c(n.categories,n.dimensions,n.dimensions))
for(j in 1:n.categories){
  sigma0[j,,]<-rWishart(1,n.dimensions*n.dimensions,Posdef(n.dimensions))
}
alpha<-rep(1,n.categories)






category.probability<-matrix(ncol=n.categries,nrow=n.run)

MCMC<-function(Y,mu0,Lambda0,nu0,sigma0,alpha){
  n.categories<-2
  n<-dim(Y)[1]
  n.dimensions<-dim(Y)[2]
  
  ybar<-numeric(n.categories)
  
  for(j in 1:n.categories){
    ybar[j]<-mean(Y[,j])
  }
  
  ##MCMC setup##
  nrun=10000
  thin=
    burn=
    
    
    
  Z<-array(dim=c(nrun,n,n.categories))
  theta<-array(dim=c(nrun,n.categories))
  Mu<-array(dim=c(nrun,n.categories,n.dimensions))
  Sigma<-array(dim=c(nrun,n.categories,n.dimensions,n.dimensions))

  ##set initial values##
  for(s in 1:n){
  Z[1,s,]<-rmultinom(1,1,rep(1/n.categories,n.categories))
  }
  theta[1,]<-rep(1/n.categories,n.categories)
  for(j in 1:n.categories){
  Mu[1,j,]<-rmvnorm(mu0,Lambda0)
  Sigma[1,j,,]<-solve(rWishart(1,nu0,sigma0))
  }
  
  ##begin sampler##
  for(iter in 2:nrun){
    
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
          ?rmultinom
          Z.jump[s,]<-rmultinom(1,1,prob=P)
        }
    loglik.current<-numeric(n)
    loglik.jump<-numeric(n)
    prior.current<-numeric(n)
    prior.jump<-numeric(n)
    
    for(s in 1:n){
      j.jump<-which(Z.jump[s,]==1)
      j.current<-which(Z[(iter-1),s,]==1)
      loglik.current[s]<-log(dmvnorm(Y[s,],Mu[iter,j.current,],Sigma[iter,j.current,,]))
      loglik.jump[s]<-log(dmvnorm(Y[s,],Mu[iter,j.jump],Sigma[iter,j.jump,,]))
      prior.current[s]<-log(dmultinom(Z[(iter-1),s,],1,theta[(iter-1),]))
      prior.jump[s]<-log(dmultinom(Z.jump[s,]),1,theta[(iter-1),])
    }
    
    accept.prob<- min(1, exp(sum(loglik.jump+prior.jump))/exp(sum(loglik.current+prior.current)))
    accept<-rbinom(1,1,accept.prob)
    if(accept == 1){
      Z[iter,,]<-Z.jump
    }
    if(accept==0){
      Z[iter,,]<-Z[(iter-1),,]
    }
    
    ##update theta##
    Z.sum<-numeric(n.categories)
    
    for(j in 1:n.categories){
      Z.sum[j]<-sum(Z[iter,,j])
    }
    
    theta[iter,]<-theta[(iter-1),]+Z.sum
    
    ##update mu##
    for(j in 1:m){
      Lambdan<-solve(solve(Lambda0)+n*solve(Sigma[(iter-1),j,,]))
      mun<-Lambdan%*%(solve(Lambda0)*mu0+n*solve(Sigma[(iter-1),j,,])%*%ybar[j])
      Mu[s,j,]<-rmvnorm(1,mun,Lambdan)
    }
    
    ##update sigma##
    for(j in 1:m){
      nun<-nu0 + n
      Sn<-Sigma0+(t(Y[,,j])-Mu[(iter-1),j,])%*%t(t(Y[,,j])-Mu[(iter-1),j,])
      Sigma[iter,j,,]<-solve(rwishart(1,nun,solve(Sn)))
    }
  }  
}
MCMC(Y,mu0,Lambda0,nu0,sigma0,alpha)
  ##weight matrix stuff##
  
  for(j in 1:n.categories){
    p<-matrix(nrow=n.dimensions,ncol=n.dimensions)
    for(d in 1:n.dimensions){
      p[d,]<-rgamma(n.dimensions,alpha[d,])
      for(i in 1:n.dimensions){                                                                                        
  
        
  
        A[s,j,d,i]<-p[d,i]/sum(p[d,])
      }
  
    }
  }
}
  
MCMC(Y,mu0,Lambda0,nu0,sigma0,alpha)
