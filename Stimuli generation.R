###This probably won't be used##

feature.generator.l.r<-function(n.pixels){
  segmant<-matrix(nrow=n.pixels,ncol=n.pixels)
  for(i in 1:n.pixels){
    for(j in 1:n.pixels){
      segmant[i,j]<-rbinom(1,1,.2)
      if(i==1 && j==1){
        segmant[i,j]<-1  
      }
      if(i==n.pixels && j==n.pixels){
        segmant[i,j]<- 1
      }  
    }
  }
  for(i in 1:(n.pixels-1)){
    for(j in 1:(n.pixels-1)){
      if(segmant[i,j]==1){
        if(segmant[i,(j+1)]==0 || segmant[(i+1),(j+1)]==0 || segmant[(i+1),j]==0){
          w<-sample.int(1, 1:3)
          if(w==1){
            segmant[i,(j+1)]<-1
          }
          if(w==2){
            segmant[(i+1),j]<-1
          }
          if(w==3){
            segmant[(i+1),(j+1)]<-1
          }
        }
      }  
    }
  }
  return(segmant)
}

feature.generator.r.l<-function(n.pixels){
  segmant<-matrix(nrow=n.pixels,ncol=n.pixels)
  for(i in 1:n.pixels){
    for(j in 1:n.pixels){
      segmant[i,j]<-rbinom(1,1,.3)
      if(i==1 && j==n.pixels){
        segmant[i,j]<-1  
      }
      if(i==n.pixels && j==1){
        segmant[i,j]<- 1
      }  
    }
  }
  for(i in 1:(n.pixels-1)){
    for(j in 2:(n.pixels)){
      if(segmant[i,j]==1){
        if(segmant[i,(j-1)]==0 || segmant[(i+1),(j-1)]==0 || segmant[(i+1),j]==0){
          w<-sample.int(1, 1:3)
          if(w==1){
            segmant[i,(j-1)]<-1
          }
          if(w==2){
            segmant[(i+1),j]<-1
          }
          if(w==3){
            segmant[(i+1),(j-1)]<-1
          }
        }
      }  
    }
  }
  return(segmant)
}

hello<-feature.generator.l.r(4)
print(hello)
image(hello,col=gray(12:1/12))


Features<-list()
Features[[8]]<-hello

create.stimuli<-function(f1,f2){
  for(o in 1:length(f1)){
    stimuli<-matrix(rep(length(f1)*length(f2)),nrow=length(f1),nrow=length(f2))
    values<-f1[[o]]
    for(i in 1:length(values[,1])){
      for(j in 1:length(values[1,])){
        stimuli[,o*i]
      }
    }
  }
}

image(matrix(m[i,], nrow=28)[,28:1], 
      col=gray(12:1/12), # color scheme for visualizing
      xlab=which(mL[i,]==1)-1, # extract the label from the label matrix
      xaxt='n',yaxt='n' # remove x and y axes.
)
?image
