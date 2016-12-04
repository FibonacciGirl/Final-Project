#install.packages("mvtnorm")
#install.packages("MCMCpack")
#instal.packages("tmvtnorm")
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

  


#latent class stuff#

MCMC<-function(Y,ybar,mu0,Lambda0,nu0,sigma0,alpha,W0){
  n.categories<-2
  n<-dim(Y)[2]
  n.dimensions<-dim(Y)[1]
  
  ##set initial values##
  Z<-array(dim=c(nrun,n,n.categories))
  theta<-array(dim=c(nrun,n.categories))
  Mu<-array(dim=c(nrun,n.categories,n.dimensions))
  Sigma<-array(dim=c(nrun,n.categories,n.dimensions,n.dimensions))
  W<-array(dim=c(nrun,n.categories,n.dimensions,n.dimensions))
  
  theta[1,]<-rep(1/n.categories,n.categories)
  for(s in 1:n){
    Z[1,s,]<-rmultinom(1,1,theta[1,])
  }
  for(j in 1:n.categories){
    Mu[1,j,]<-rmvnorm(1,mu0[j,],Lambda0[j,,])
    W[1,j,,]<-rWishart(1,n.dimensions,W0[j,,])
    Sigma.inv<-rWishart(1,n.dimensions,sigma0[j,,])
    Sigma[1,j,,]<-solve(matrix(Sigma.inv,nrow=n.dimensions,ncol=n.dimensions))
  }
  
  
  ##MCMC setup##
  nrun =10000
  burn = 5000
  thin = 10
  eff.sam = (nrun-burn)/thin
  
  Theta.out<-array(dim=c(eff.sam,n.categories))
  Sigma.out<-array(dim=c(eff.sam,n.categories,n.dimensions,n.dimensions))
  Mu.out<-array(dim=c(eff.sam,n.categories,n.dimensions))
  Z.out<-array(dim=c(eff.sam,n,n.categories))
  W.out<-array(dim=c(eff.sam,n.categories,n.dimensions,n.dimensions))
  
  ##begin sampler##
  for(iter in 2:nrun){
    print(iter)
    ##update Z and W##
    
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
    
    W.jump<-array(dim=c(n.categories,n.dimensions,n.dimensions))
    
    for(j in 1:n.categories){
      W.jump[j,,]<-rWishart(1,n.dimensions,W0[j,,])
    }
    
    loglik.w.current<-numeric(n)
    loglik.w.jump<-numeric(n)
    prior.w.current<-numeric(n)
    prior.w.jump<-numeric(n)
    
    loglik.z.current<-numeric(n)
    loglik.z.jump<-numeric(n)
    prior.z.current<-numeric(n)
    prior.z.jump<-numeric(n)
    
    for(s in 1:n){
      j.z.jump<-which(Z.jump[s,]==1)
      j.z.current<-which(Z[(iter-1),s,]==1)
      loglik.z.current[s]<-log(dmvnorm(Y[,s],W[(iter-1),j.z.current,,]%*%Mu[(iter-1),j.z.current,],W[(iter-1),j.z.current,,]%*%Sigma[(iter-1),j.z.current,,]%*%t(W[(iter-1),j.z.current,,])))
      loglik.z.jump[s]<-log(dmvnorm(Y[,s],W[(iter-1),j.z.jump,,]%*%Mu[(iter-1),j.z.jump,],W[(iter-1),j.z.jump,,]%*%Sigma[(iter-1),j.z.jump,,]%*%t(W[(iter-1),j.z.jump,,])))
      prior.z.current[s]<-log(dmultinom(Z[(iter-1),s,],1,theta[(iter-1),]))
      prior.z.jump[s]<-log(dmultinom(Z.jump[s,],1,theta[(iter-1),]))
    
      
      j.w.jump<-which(Z[(iter-1),s,]==1)
      j.w.current<-which(Z[(iter-1),s,]==1)
      ##get rid of repeat
      
      loglik.w.current[s]<-log(dmvnorm(Y[,s],W[(iter-1),j.w.current,,]%*%Mu[(iter-1),j.w.current,],W[(iter-1),j.w.current,,]%*%Sigma[(iter-1),j.w.current,,]%*%t(W[(iter-1),j.w.current,,])))
      loglik.w.jump[s]<-log(dmvnorm(Y[,s],W.jump[j.w.jump,,]%*%Mu[(iter-1),j.w.jump,],W.jump[j.w.jump,,]%*%Sigma[(iter-1),j.w.jump,,]%*%t(W.jump[j.w.jump,,])))
      prior.w.current[s]<-log(dwish(W[(iter-1),j.w.current,,],n.dimensions,W0[j,,]))
      prior.w.jump[s]<-log(dwish(W.jump[j.w.jump,,],n.dimensions,W0[j,,]))
      }
    
    
    prob.z.jump<-sum(loglik.z.jump)+sum(prior.z.jump)
    prob.z.current<-sum(loglik.z.current)+sum(prior.z.current)
    
    accept.z.prob<- min(1,exp(prob.z.jump)/max(.00001,exp(prob.z.current)))
    
    accept.z<-rbinom(1,1,accept.z.prob)
    
    if(accept.z == 1){
      Z[iter,,]<-Z.jump
    }
    if(accept.z == 0){
      Z[iter,,]<-Z[(iter-1),,]
    }
    
    prob.w.jump<-sum(loglik.w.jump)+sum(prior.w.jump)
    prob.w.current<-sum(loglik.w.current)+sum(prior.w.current)
    
    accept.w.prob<- min(1,exp(prob.w.jump)/max(.00001,exp(prob.w.current)))
    
    accept.w<-rbinom(1,1,accept.w.prob)
    
    if(accept.w == 1){
      W[iter,,,]<-W.jump
    }
    if(accept.w == 0){
      W[iter,,,]<-W[(iter-1),,,]
    }
    ##update theta##
    Z.sum<-numeric(n.categories)
    alpha_n<-numeric(n.categories)
    
    for(j in 1:n.categories){
      Z.sum[j]<-sum(Z[(iter-1),,j])
      alpha_n[j]<-alpha[j] +Z.sum[j]
    }
    theta[iter,]<-rdirichlet(1,alpha_n)
    
    #print(theta[iter,])
    
    ##update mu##
    tau<-array(dim=c(n.categories,n.dimensions,n.dimensions))
    for(j in 1:n.categories){
      tau[j,,]<-solve(W[(iter-1),j,,]%*%Sigma[(iter-1),j,,]%*%t(W[(iter-1),j,,]))
      print(tau[j,,])
      
      A.inv<-nj[j]*(t(W[(iter-1),j,,])%*%tau[j,,]%*%W[iter-1,j,,]) + solve(Lambda0[j,,])

      b<-solve(Lambda0[j,,])%*%mu0[j,] + nj[j]*t(W[(iter-1),j,,])%*%tau[j,,]%*%ybar[j,]
      Mu[iter,j,]<-rmvnorm(1,A.inv%*%b,A.inv)
    }
    
    ##update sigma##
  
        
      nun<-n.dimensions + n
    for(j in 1:n.categories){ 
      Y.j<-array(dim=c(n.categories,n.dimensions,nj[j]))
      Y.j[j,,]<-Y[,which(category==j)]
      Y.w.j<-array(dim=c(nj[j],n.dimensions,n.dimensions))
      Sig<-array(dim=c(nj[j],n.dimensions,n.dimensions))
      Sw<-array(dim=c(n.categories,n.dimensions,n.dimensions))
      for(i in 1:nj[j]){
        Y.w.j[i,,]<-(Y.j[j,,i]-W[(iter-1),j,,]%*%Mu[(iter-1),j,])%*%t(Y.j[j,,i]-W[(iter-1),j,,]%*%Mu[(iter-1),j,])
        Sig[i,,]<- solve(W[(iter-1),j,,])%*%(Y.w.j[i,,])%*%solve(t(W[(iter-1),j,,]))    
      }
      Sw[j,1,1]<-sum(Sig[,1,1])
      Sw[j,1,2]<-sum(Sig[,1,2])
      Sw[j,2,1]<-sum(Sig[,2,1])
      Sw[j,2,2]<-sum(Sig[,2,2])

      sig.inv<-rWishart(1,nun,solve(Sw[j,,]))
      Sigma[iter,j,,]<-solve(matrix(sig.inv,nrow=n.dimensions,ncol=n.dimensions))
    }
    
    
    ##burn and thin##
    if(iter %% thin == 0 && iter>burn){
      sav=(iter-burn)/thin
      Theta.out[sav,]<-theta[iter,]
      Sigma.out[sav,,,]<-Sigma[iter,,,]
      Mu.out[sav,,]<-Mu[iter,,]
      Z.out[sav,,]<-Z[iter,,]
      W.out[sav,,,]<-W[iter,,,]
    }
  }
  #return(theta,Sigma, Mu,Z)
  return(list(THETA=Theta.out,SIGMA=Sigma.out,MU=Mu.out,Z.OUT=Z.out,W.OUT=W.out))
}




##data##
data<-read.csv('TD')
category<-data$category
data$category<-NULL
n.categories<-length(unique(category))
Y<-t(data.matrix(data))


n<-dim(Y)[2]
n.dimensions<-dim(Y)[1]

ybar<-array(dim= c(n.categories,n.dimensions))
nj<-numeric()

for(j in 1:n.categories){
  ybar[j,]<-mean(Y[,which(category==j)])
  nj[j]<-length(Y[,which(category==j)])
}
#choose fixed hyper paratments#
mu0<-array(dim=c(n.categories,n.dimensions))
mu0<-rmvnorm(n.categories,rep(0,n.dimensions),Posdef(n.dimensions))
Lambda0<-array(dim=c(n.categories,n.dimensions,n.dimensions))
for(j in 1:n.categories){
  Lambda0[j,,]<-rWishart(1,n.dimensions,Posdef(n.dimensions))
}
Lambda0[1,,]

nu0<-numeric(n.categories)
nu0<-rep(n.dimensions,n.categories)
sigma0<-array(dim=c(n.categories,n.dimensions,n.dimensions))
for(j in 1:n.categories){
  sigma0[j,,]<-rWishart(1,n.dimensions,Posdef(n.dimensions))
}
W0<-array(dim=c(n.categories,n.dimensions,n.dimensions))
for(j in 1:n.categories){
  W0[j,,]<-rWishart(1,n.dimensions,Posdef(n.dimensions))
}
alpha<-rep(1,n.categories)

test<- MCMC(Y,ybar,mu0,Lambda0,nu0,sigma0,alpha,W0)


th<-test[[1]]
si<-test[[2]]
mu<-test[[3]]
z.z<-test[[4]]
w.w<-test[[5]]

density(mu)


plot(density(w.w[,2,2,1]))
plot(density(z.z[,1,]))
sigma<-as.mcmc.list(sigma1.mcmc,sigma2.mcmc)
traceplot(sigma1.mcmc)
