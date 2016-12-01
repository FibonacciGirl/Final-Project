c<-10 #slope of exp decay function
concept.bias<-1
learning.rate.catergory.weights<-.01 #for updating category specific weights
learning.rate.exemplar.weight<-.1 #for updating exemplar specific weights
weight.combination<-.5 #combining specific weights



similarity.weights<-function(td,stim){

  exemplars<-data.matrix(td)
  n.dimensions<-dim(exemplars)[1]
  n.exemplars<-dim(exemplars)[2]
  print(n.exemplars)
  distance<-matrix(nrow=n.exemplars,ncol=n.dimensions)
  weights<-matrix(nrow=n.exemplars,ncol=n.dimensions)
  normal.weights<-matrix(nrow=n.exemplars,ncol=n.dimensions)
  index<-matrix(nrow=n.exemplars,ncol=n.dimensions)


  #for each exemplar create distance vector that store distance to the stimuli on dimension d
  
  for(e in 1:n.exemplars){ 
    #for each dimension
    for (d in 1:n.dimensions){
      distance[e,d]<-sqrt((exemplars[e,d]-stim[d])^2)
    }
  #sort distance in decreasing order and also return the index of the exemplars associated with the ranks 
    weight<-sort(distance[,d],decreasing=TRUE,index.return=TRUE)
    index[,d]<-weight[[2]]
    weights[,d]<-weight[[1]]
  #change weights to reflect rank, including ties, and adding noise
    mu<-10
      for(m in 1:max(weights[,d])){
        tie<-1
        if(any(duplicated(weights[,d]))){
          tie<-length(which(weights[,d]==m)) 
        }
        for(e in 1:n.exemplars){
          if(weights[e,d]==m){
            sd<-sd*tie
            weights[e,d]<-rnorm(1,mu,sd)
            while(weights[e,d]<0){
              weights[e,d]<-rnorm(1,mu,sd)
            }
          }
        }
      mu<-mu+1
      }
    #normalize across dimension
    total<-sum(weights[,d])
    for(e in 1:n.exemplars){
      weights[e,d]<-weights[e,d]/total
    }
  }  
  #normalize so the sum of the dimension weights in a particular exemplar eqauls 1
  for(e in 1:n.exemplars){
    for(d in 1:n.dimensions){
      normal.weights[e,d]<-weights[e,d]/sum(weights[e,])
    }
  }
  return(cbind(normal.weights,index))
}

similarity<-function(distance){
  c<-1
  similarity<-exp(-c*distance)
}

probability<-function(td,stim,cat){
  category<-c(cat)
  exemplars<-data.matrix(td)
  n.categories<-length(unique(cat))
  n.dimensions<-dim(exemplars)[2]
  n.exemplars<-dim(exemplars)[1]  
  w<-similarity.weights(td,stim) #create weights
  
  similar<-numeric(n.categories) #to hold similarity for each category
  prob<-numeric(n.categories) #to hold probability for each category

  #go through each of the categories and get the similarity value for the stimuli
  for(c in 1:n.categories){
    distance<-c()
    for(e in 1:n.exemplars){
      if(category[e] == c){  
        for(d in 1:n.dimensions){
   
        distance<-c(distance,sqrt(w[which(w[,n.dimensions+d]==e),d]*(exemplars[e,d]-stim[d])^2))
        #select proper weight by the index stored by weighting function
        }
      }
    }
    print(as.vector(distance))
    similar[c]<-sum(sapply(distance,similarity))
  }
  
  #conver similarity values to probability
  total.similar<-sum(similar)
  for(c in 1:n.categories){
    prob[c]<-similar[c]/total.similar
  }
  print(prob)
  return(prob)
}

