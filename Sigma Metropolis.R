install.packages('MCMCpack')
library('MCMCpack')

?dwish

Sigma.jump<-array(dim=(n.categories,n.dimensions,n.dimensions))
for(j in 1:n.categories){
  Sigma.jump[j,,]<-rWishart(1,n.dimensions,W0)
}

loglik.Sigma.current<-numeric(n)
loglik.Sigma.jump<-numeric(n)
prior.Sigma.current<-numeric(n)
prior.Sigma.jump<-numeric(n)

for(j in 1:n){
  j.Sigma.jump<-which(Z[s,]==1)
  j.Sigma.current<-which(Z[(iter-1),s,]==1)
  for(s in 1:n){
    loglik.Sigma.current[s]<-log(dmvnorm(Y[,s],W[(iter-1),j.Sigma.current,,]%*%Mu[(iter-1),j.Sigma.current,],W[(iter-1),j.Sigma.current,,]%*%Sigma[(iter-1),j.Sigma.current,,]%*%t(W[(iter-1),j.Sigma.current,,])))
    loglik.Sigma.jump[s]<-log(dmvnorm(Y[,s],W[(iter-1),j.Sigma.jump,,]%*%Mu[(iter-1),j.Sigma.jump,],W[(iter-1),j.Sigma.jump,,]%*%Sigma[(j.Sigma.jump,,]%*%t(W[(iter-1),j.Sigma.jump,,])))
    prior.Sigma.current[s]<-log(dwish(solve(Sigma[(iter-1),j.Sigma.current,,]),n.dimensions,S0))
    prior.Sigma.jump[s]<-log(dwish(Sigma.jump[j.Sigma.jump,,],n.dimensions,S0))
  }
}

prob.Sigma.jump<-sum(loglik.Sigma.jump)+sum(prior.Sigma.jump)
prob.Sigma.current<-sum(loglik.Sigma.current)+sum(prior.Sigma.current)

accept.Sigma.prob<- min(1,exp(prob.Sigma.jump)/max(.00001,exp(prob.Sigma.current)))

accept.Sigma<-rbinom(1,1,accept.Sigma.prob)

if(accept.Sigma == 1){
  Sigma[iter,,,]<-Sigma.jump
}
if(accept.Sigma == 0){
  Sigma[iter,,,]<-Sigma[(iter-1),,,]
}