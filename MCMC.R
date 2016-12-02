
n.categories<-4 # run the loop this many times to get prob for each category
n.dimensions<-3 #number of dimensions in which the stimulus is observed
training.data<-
test.data<-


  
  #latent class stuff#

##data##
  
Y<-array(n,n.dimensions)
ybar<-numeric(n.categories)

  
for(j in 1:n.categories){
  ybar[j]<-mean(Y[,,j])
}

##MCMC setup##
nrun=10000
thin=
burn=
  
  

Z<-array(dim=c(nrun,n.categories))
Mu<-array(dim=c(nrun,n.categories,n.dimensions))
Sigma<-array(dim=c(nrun,n.categories,n.dimensions,n.dimesnsions))
A<-array(dim=c(nrun,n.categories,n.dimensions,n.dimensions))

mu0<-array(dim=c(n.categories,n.dimensions))
Lambda0<-array(dim=c(n.categories,n.dimensions,n.dimensions))
nu0<-array(dim=c(n.categories,n.dimensions))
Sigma0<-array(dim=c(n.categories,n.dimensions))

alpha<-array(dim=c(n.categories))  


##latent class##
Z<-array(dim=c(nrun,n,n.categories))

alpha.z<-array(dim=c(n.run,n.categories))


category.probability<-matrix(ncol=n.categries,nrow=n.run)

Z.jump<-array(dim=c(n,n.categories))

    for(s in 1:n){
      X<-numeric(n.categories)
      for(i in 1:n.categories){
        X[i]<-rgamma(1,alpha[(iter-1),i])
      }
      X.sum<-sum(X)
      P<-numeric(n.categories)
      for(i in 1:n.categories){
        P[i]<-X[i]/X.sum
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
  prior.current[s]<-log(dmultinom(Z[(iter-1),s,],1,alpha[(iter-1),]))
  prior.jump[s]<-log(dmultinom(Z.jump[s,]),1,alpha[(iter-1),])
}

accept.prob<- min(1, exp(sum(loglik.jump+prior.jump))/exp(sum(loglik.current+prior.current)))
accept<-rbinom(1,1,accept.prob)
if(accept == 1){
  Z[iter,,]<-Z.jump
}
if(accept==0){
  Z[iter,,]<-Z[(iter-1),,]
}

Z.sum<-numeric(n.categories)

for(j in 1:n.categories){
  Z.sum[j]<-sum(Z[iter,,j])
}

alpha[iter,]<-alpha[(iter-1),]+Z.sum


for(j in 1:m){
  Lambdan<-solve(solve(Lambda0)+n*solve(Sigma[(s-1),j,,]))
  mun<-Lambdan%*%(solve(Lambda0)*mu0+n*solve(Sigma[(s-1),j,,])%*%ybar[j])
  Mu[s,j,]<-rmvnorm(1,mun,Lambdan)
}

for(j in 1:m){
  nun<-nu0 + n
  Sn<-Sigma0+(t(Y[,,j])-Mu[(s-1),j,])%*%t(t(Y[,,j])-Mu[(s-1),j,])
  Sigma[s,j,,]<-solve(rwishart(1,nun,solve(Sn)))
}

for(j in 1:n.categories){
  p<-matrix(nrow=n.dimensions,ncol=n.dimensions)
  for(d in 1:n.dimensions){
    p[d,]<-rgamma(n.dimensions,alpha[(s-1),d,])
    for(i in 1:n.dimensions){                                                                                        

      
##weight matrix stuff##
      A[s,j,d,i]<-p[d,i]/sum(p[d,])
    }

  }
}