install.packages('MCMCpack')
library('MCMCpack')

?dwish
##full conditionals##
for(j in 1:n.categories){
A.inv<-n*(t(W[(iter-1),j,,]%*%solve(W[(iter-1),j,,]%*%Sigma[(iter-1),j,,]%*%t(W[(iter-1),j,,])))) + solve(Lambda0[j,,])
b<-solve(Lambda0)%*%mu0[j,] + n*solve(W[(iter-1),j,,]%*%Sigma[(iter-1),j,,]%*%t(W[(iter-1),j,,]))%*%ybar[j,]
Mu[iter,j,]<-rmvnorm(1,A%*%b,A)
}

Y.w.j<-numeric(n.categories)
for(j in 1:n.categories){
  nun<-n.dimensions + n
  Y.w.j[j]<-sum((Y[,which(category==j)]-W[(iter-1),j,,]%*%Mu[(iter-1),j,,])%*%t(Y[,which(category==j)]-W[(iter-1),j,,]%*%Mu[(iter-1),j,,]))
  
  Sw<-solve(W[(iter-1),j,,])%*%(Y.w.j[j]*solve(t(W[(iter-1),j,,])))
  
  sig.inv<-rWishart(1,nun,solve(Sw))
  Sigma[iter,j,,]<-solve(matrix(sig.inv,nrow=n.dimensions,ncol=n.dimensions))
}


##metropolis##

W.jump<-array(dim=(n.categories,n.dimensions,n.dimensions))
for(j in 1:n.categories){
  W.jump[j,,]<-rWishart(1,n.dimensions,W0)
}

loglik.w.current<-numeric(n)
loglik.w.jump<-numeric(n)
prior.w.current<-numeric(n)
prior.w.jump<-numeric(n)

for(s in 1:n){
  j.w.jump<-which(Z[s,]==1)
  j.w.current<-which(Z[(iter-1),s,]==1)
    loglik.w.current[s]<-log(dmvnorm(Y[,s],W[(iter-1),j.w.current,,]%*%Mu[(iter-1),j.w.current,],W[(iter-1),j.w.current,,]%*%Sigma[(iter-1),j.w.current,,]%*%t(W[(iter-1),j.w.current,,])))
    loglik.w.jump[s]<-log(dmvnorm(Y[,s],W.jump[j.w.jump,,]%*%Mu[(iter-1),j.w.jump,],W.jump[j.w.jump,,]%*%Sigma[(iter-1),j.w.jump,,]%*%t(W.jump[j.w.jump,,])))
    prior.w.current[s]<-log(dwish(W[(iter-1),j.w.current,,],n.dimensions,W0))
    prior.w.jump[s]<-log(dwish(W.jump[j.w.jump,,],n.dimensions,W0))
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

