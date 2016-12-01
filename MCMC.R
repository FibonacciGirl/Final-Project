
n.categories<-4 # run the loop this many times to get prob for each category
n.dimensions<-3 #number of dimensions in which the stimulus is observed
training.data<-
test.data<-


  
  #latent class stuff#
Z<-array(dim=c(n.run,n,n.categories))
alpha.z<-array(dim=c(n.run,n,n.categories))

mu<-matrix(ncol=n.dimensions,nrow=n.run)
Sigma<-array(dim=c(n.dimensions,n.dimensions,n.run))

category.probability<-matrix(ncol=n.categries,nrow=n.run)

for(iter in 1:n.run){
  for(s in 1:n){
    Y<-c()
    for(i in 1:n.categories){
      Y[i]<-rgamma(1,alpha[(iter-1),i])
    }
    Y.sum<-sum(Y)
    P<-numeric(n.categories)
    for(i in 1:n.categories){
      P[i]<-Y[i]/Y.sum
    }
  Z[iter,s,]<-rmultinom(1,n.categories,P)
  alpha[iter,s,]<-alpha.z((iter-1),s,)+sum(Z[iter,s,])}
  }
}



##MCMC setup##
nrun=10000
thin=
burn=
  
Y<-#data

Z<-array(dim=c(nrun,n.categories))
Mu<-array(dim=c(nrun,n.categories,n.dimensions))
Sigma<-array(dim=c(nrun,n.categories,n.dimensions,n.dimesnsions))
A<-array(dim=c(nrun,n.categories,n.dimensions,n.dimensions))

mu0<-array(dim=c(n.categories,n.dimensions))
Lambda0<-array(dim=c(n.categories,n.dimensions,n.dimensions))
nu0<-array(dim=c(n.categories,n.dimensions))
Sigma0<-array(dim=c(n.categories,n.dimensions))

alpha<-array(dim=c(n.categories))  

for(j in 1:m){
  ybar[j]<-mean(Y[,,j])
}

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
      A[s,j,d,i]<-p[d,i]/sum(p[d,])
    }

  }
}