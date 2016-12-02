c<-1 #slope of exp decay function
concept.bias<-1
learning.rate.catergory.weights<-.01 #for updating category specific weights
learning.rate.exemplar.weight<-.1 #for updating exemplar specific weights
weight.combination<-.5 #combining specific weights



similarity.weights<-function(td,stim){

  exemplars<-data.matrix(td)
  n.dimensions<-dim(exemplars)[2]
  n.exemplars<-dim(exemplars)[1]
  distance<-matrix(nrow=n.exemplars,ncol=n.dimensions)
  weights<-matrix(nrow=n.exemplars,ncol=n.dimensions)
  normal.weights<-matrix(nrow=n.exemplars,ncol=n.dimensions)
  index<-matrix(nrow=n.exemplars,ncol=n.dimensions)
  sd<-20

  #for each dimension
  for (d in 1:n.dimensions){
    
  #for each exemplar create distance vector that store distance to the stimuli on dimension d
  for(e in 1:n.exemplars){ 
    distance[e,d]<-sqrt((exemplars[e,d]-stim[d])^2)
  }
  #sort distance in increasing order and also return the index of the exemplars associated with the ranks 
    weight<-sort(distance[,d],decreasing=FALSE,index.return=TRUE)
    index[,d]<-weight[[2]]
    weights[,d]<-weight[[1]]
  #change weights to reflect rank, including ties, and adding noise
    mu<-10
      for(m in 1:n.exemplars){
        tie<-1
        if(any(duplicated(weights[,d]))){
          tie<-length(which(weights[,d]== weights[m,d]))
        if(tie==0){
            tie<-1
          }
        }
        sd<-sd*tie
        weights[m,d]<-rnorm(1,mu,sd)
        while(weights[e,d]<0){
          weights[m,d]<-rnorm(1,mu,sd)
        }
        if(length(duplicated(weights[m,d]==0))){
          mu<-mu+10
        }
      }
    #normalize across dimension
    total<-sum(weights[,d])
    for(e in 1:n.exemplars){
      weights[e,d]<-weights[e,d]/total
    }
  }  
  #normalize so the sum of the dimension weights in a particular exemplar eqauls 1
  for(e in 1:n.exemplars){
    total<-sum(weights[e,])
    for(d in 1:n.dimensions){
      normal.weights[e,d]<-weights[e,d]/total
    }
  }
  print(sum(normal.weights))
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
  
  print (w)
    
  similar<-numeric(n.categories) #to hold similarity for each category
  prob<-numeric(n.categories) #to hold probability for each category

  #go through each of the categories and get the similarity value for the stimuli
  for(c in 1:n.categories){
    distance<-numeric(n.exemplars)
    dim<-c()
    for(e in 1:n.exemplars){
      if(category[e] == c){  
        for(d in 1:n.dimensions){
        dim<-c(dim,sqrt(w[which(w[,n.dimensions+d]==e),]*(exemplars[e,d]-stim[d])^2))
        #select proper weight by the index stored by weighting function
        }
        distance[e]<-sum(sapply(dim,similarity))  
      }
    }
    distance<-distance[!is.na(distance)]
    similar[c]<-sum(distance)
  }
  #conver similarity values to probability
  total.similar<-sum(similar)
  for(c in 1:n.categories){
    prob[c]<-similar[c]/total.similar
  }
  return(prob)
}

