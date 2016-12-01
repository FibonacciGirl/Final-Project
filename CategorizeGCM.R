GCM<-function(td){
  #convert input
  td<-data.frame(td)
  category<-as.vector(td$category)
  td$category<-NULL
  exemplars<-data.matrix(td)
  n.categories<-length(unique(category))
  n.dimensions<-dim(exemplars)[2]
  n.exemplars<-dim(exemplars)[1]
  print(n.exemplars)
  category.guess<-numeric(n.exemplars)
  for(e in 3:n.exemplars){
    td<-exemplars[(1:(e-1)),]
    stim<-c(exemplars[e,])
    cat<-category[1:(e-1)]
    prob.categorize<-probability(td,stim,cat)
    result<-rmultinom(1,1,prob.categorize)
    category.guess[e]<-which(result==1)
  }
  td<-data.frame(exemplars)
  td$category<-category
  td$category.guess<-category.guess
  return(td)
}

GCM(training.data)
