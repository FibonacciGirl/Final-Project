
n.categories<-4 # run the loop this many times to get prob for each category
n.dimensions<-3 #number of dimensions in which the stimulus is observed
training.data<-
test.data<-


  
  #latent class stuff#




##MCMC setup##
nrun=10000
thin=
burn=
  
  
Y<-array(n,n.dimensions)

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
      Y<-c()
      for(i in 1:n.categories){
        Y[i]<-rgamma(1,alpha[(iter-1),i])
      }
      Y.sum<-sum(Y)
      P<-numeric(n.categories)
      for(i in 1:n.categories){
        P[i]<-Y[i]/Y.sum
      }
      ?rmultinom
      Z.jump[s,]<-rmultinom(1,1,prob=P)
    }
loglik.current<-numeric(n.categories)
loglik.jump<-numeric(n.categories)
for(s in 1:n){
  j.jump<-which(Z.jump[s,]==1)
  j.current<-which(Z[(iter-1),s,]==1)
  loglik.current[s]->log(dmvnorm(Y[s,],Mu[iter,j.current,],Sigma[iter,j.current,,]))
  loglik.jump[s]->(dmvnorm(Y[s,],Mu[iter,j.jump],Sigma[iter,j.jump,,]))
}

accept.prob<- min(1, exp(sum(loglik.jump))/exp(sum(loglik.current)))
accept<-rbinom(1,1,accept.prob)
if(accept == 1){
  Z[iter,,]<-Z.jump
}
if(accept==0){
  Z[iter,,]<-Z[(iter-1),,]
}

alpha[s,]<-alpha.z[(iter-1),s,]+sum(Z[iter,s,])

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

      
##weight matrix stuff##
      A[s,j,d,i]<-p[d,i]/sum(p[d,])
    }

  }
}