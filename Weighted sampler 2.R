#install.packages("mvtnorm")
#install.packages("MCMCpack")
#instal.packages("tmvtnorm")
#install.packages("Matrix")
library(Matrix)
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

MCMC<-function(Y){
  
  n.categories<-length(unique(category))
  n<-dim(Y)[2]
  n.dimensions<-dim(Y)[1]
  
  ybar<-array(dim= c(n.categories,n.dimensions))
  nj<-numeric()
  
  for(j in 1:n.categories){
    ybar[j,]<-mean(Y[,which(category==j)])
    nj[j]<-length(Y[,which(category==j)])
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
  
  
  #choose fixed hyper paratments#
  mu0<-array(rep(0,n.dimensions*n.categories),dim=c(n.categories,n.dimensions))
  Lambda0<-array(dim=c(n.categories,n.dimensions,n.dimensions))
  for(j in 1:n.categories){
    Lambda0[j,,]<-rWishart(1,n.dimensions,diag(1,n.dimensions,n.dimensions))
  }
  Lambda0[1,,]
  
  
  nu0<-rep(n.dimensions,n.categories)
  sigma0<-rWishart(1,n.dimensions,diag(100,n.dimensions,n.dimensions))
  sigma0<-solve(matrix(sigma0,nrow=n.dimensions,ncol=n.dimensions))
  
  W0<-array(dim=c(n.categories,n.dimensions,n.dimensions))
  for(j in 1:n.categories){
    W0[j,,]<-rWishart(1,n.dimensions,diag(1,n.dimensions,n.dimensions))
  }
  alpha<-rep(1,n.categories)
  
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
    Sigma.inv<-rWishart(1,n.dimensions,sigma0)
    Sigma[1,j,,]<-solve(matrix(Sigma.inv,nrow=n.dimensions,ncol=n.dimensions))
  }
  
  ##begin sampler##
  for(iter in 2:nrun){
    print(iter)
    witer<-array(dim=c(n.categories,n.dimensions,n.dimensions))
    wit<-array(dim=c(n.categories,n.dimensions*n.dimensions))
    for(j in 1:n.categories){
    wit[j,]<-as.vector(nearPD(W[(iter-1),j,,])$mat) 
    witer[j,,]<-matrix(wit[j,],nrow=n.dimensions,ncol=n.dimensions)
    }
    ##update Z and W##
    
    Z.jump<-array(dim=c(n,n.categories))
    for(s in 1:n){
      Z.jump[s,]<-rmultinom(1,1,prob=theta[iter-1,])
    }
    W.jump<-array(dim=c(n.categories,n.dimensions,n.dimensions))
    
    for(j in 1:n.categories){
      W.jump[j,,]<-rWishart(1,n.dimensions,witer[j,,])
    }
    
    loglik.w.current<-numeric()
    loglik.w.jump<-numeric()
    prior.w.current<-numeric()
    prior.w.jump<-numeric()
    
    loglik.z.current<-numeric()
    loglik.z.jump<-numeric()
    prior.z.current<-numeric()
    prior.z.jump<-numeric()
    
    for(s in 1:n){
      j.z.jump<-which(Z.jump[s,]==1)
      j.z.current<-which(Z[(iter-1),s,]==1)
      WW.c<-as.vector(nearPD(witer[j.z.current,,]%*%Sigma[(iter-1),j.z.current,,]%*%t(witer[j.z.current,,]))$mat)
      WW.c<-matrix(WW.c,nrow= n.dimensions,ncol=n.dimensions)
      WW.j<-as.vector(nearPD(witer[j.z.jump,,]%*%Sigma[(iter-1),j.z.jump,,]%*%t(witer[j.z.jump,,]))$mat)
      WW.j<-matrix(WW.j,nrow= n.dimensions,ncol=n.dimensions)         
                                                               
      loglik.z.current<-loglik.z.current + log(dmvnorm(Y[,s],witer[j.z.current,,]%*%Mu[(iter-1),j.z.current,],WW.c))
      loglik.z.jump<-loglik.z.jump+log(dmvnorm(Y[,s],witer[j.z.jump,,]%*%Mu[(iter-1),j.z.jump,],WW.j))
      prior.z.current<-prior.z.current+log(dmultinom(Z[(iter-1),s,],1,theta[(iter-1),]))
      prior.z.jump<-prior.z.jump+log(dmultinom(Z.jump[s,],1,theta[(iter-1),]))
        
  
      ##get rid of repeat
      jump.w.sd<-as.vector(nearPD(W.jump[j.z.jump,,]%*%Sigma[(iter-1),j.z.current,,]%*%t(W.jump[j.z.current,,]))$mat)
      jump.w.sd<-matrix(jump.w.sd,nrow=n.dimensions,ncol=n.dimensions)                  
      loglik.w.current<-loglik.w.current+log(dmvnorm(Y[,s],witer[j.z.current,,]%*%Mu[(iter-1),j.z.current,],WW.c))
      loglik.w.jump<-loglik.w.jump +log(dmvnorm(Y[,s],W.jump[j.z.current,,]%*%Mu[(iter-1),j.z.current,],jump.w.sd))
      prior.w.current<-prior.w.current+log(dwish(witer[j.z.current,,],n.dimensions,W0[j,,]))
      prior.w.jump<-prior.w.jump+log(dwish(W.jump[j.z.current,,],n.dimensions,W0[j,,]))

    }
    
    
    prob.z.jump<-loglik.z.jump+prior.z.jump
    prob.z.current<-loglik.z.current+prior.z.current
    
    accept.z.prob<- min(1,exp(prob.z.jump)/max(.00001,exp(prob.z.current)))

    accept.z<-rbinom(1,1,accept.z.prob)
    
    if(accept.z == 1){
      Z[iter,,]<-Z.jump
    }
    if(accept.z == 0){
      Z[iter,,]<-Z[(iter-1),,]
    }
    #print('############# UPDATED Z ##############')
      
    prob.w.jump = loglik.w.jump + prior.w.jump
    prob.w.current = loglik.w.current+prior.w.current
    
    accept.w.prob<- min(1,exp(prob.w.jump)/max(.00001,exp(prob.w.current)))
    
    accept.w<-rbinom(1,1,accept.w.prob)
    
    if(accept.w == 1){
      W[iter,,,]<-W.jump
    }
    if(accept.w == 0){
      W[iter,,,]<-witer
    }

    #print('#############  UPDATED W  #############')
    ##update theta##
    Z.sum<-numeric(n.categories)
    alpha_n<-numeric(n.categories)
    z.j<-numeric(n)
    for(s in 1:n){
    z.j[s]<-which(Z[(iter-1),s,]==1)
    }
    for(j in 1:n.categories){
      Z.sum[j]<-length(which(z.j==j))
      alpha_n[j]<-alpha[j] +Z.sum[j]
    }
    theta[iter,]<-rdirichlet(1,alpha_n)
    
    #print(theta[iter,])
    
    ##update mu##
    tau<-array(dim=c(n.categories,n.dimensions,n.dimensions))
    for(j in 1:n.categories){
      h<-as.vector(nearPD(witer[j,,]%*%Sigma[(iter-1),j,,]%*%t(witer[j,,]))$mat)
      h<-matrix(h,nrow=n.dimensions,ncol=n.dimensions)
      A<-nj[j]*(t(witer[j,,])%*%solve(h)%*%witer[j,,]) + solve(Lambda0[j,,])
      A<-as.vector(nearPD(A)$mat)
      A<-solve(matrix(A,nrow=n.dimensions,ncol=n.dimensions))

      b<-solve(Lambda0[j,,])%*%mu0[j,] + nj[j]*t(witer[j,,])%*%solve(h)%*%ybar[j,]
      Mu[iter,j,]<-rmvnorm(1,A%*%b,A)
    }

    #print('################   Mu updated   ################')
    ##update sigma##
    
    
    nun<-n.dimensions + n
    Sw<-array(dim=c(n.categories,n.dimensions,n.dimensions))
  for(j in 1:n.categories){
    sum.sig<-array(sigma0,dim=c(n.dimensions,n.dimensions))
    for(i in 1:n){
      if(category[i]==j){
        sum.sig<- sum.sig + solve(witer[j,,])%*%(Y[,i]-witer[j,,]%*%Mu[(iter-1),j,])%*%t(Y[,i]-witer[j,,]%*%Mu[(iter-1),j,])%*%solve(t(witer[j,,]))
      }
    }

      Sw[j,,]<-as.vector(nearPD(sum.sig)$mat)
      Sw[j,,]<-solve(matrix(Sw[j,,],nrow=n.dimensions,ncol=n.dimensions))


      sig.inv<-rWishart(1,nun,Sw[j,,])
      print(sig.inv)
      Sigma[iter,j,,]<-solve(matrix(sig.inv,nrow=n.dimensions,ncol=n.dimensions))

  }
    print(Sigma[iter,,,])

    #print('################   Sigma updated   ################')
    
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
data<-read.csv('TD10',header=TRUE)
category<-data$category
data$category<-NULL
x<-data$x
y<-data$y

Y<-t(cbind(as.matrix(x),as.matrix(y)))

n.categories<-length(unique(category))
n<-dim(Y)[2]
n.dimensions<-dim(Y)[1]

ybar<-array(dim= c(n.categories,n.dimensions))
nj<-numeric()

for(j in 1:n.categories){
  ybar[j,]<-mean(Y[,which(category==j)])
  nj[j]<-length(Y[,which(category==j)])
}


test<- MCMC(Y)


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
