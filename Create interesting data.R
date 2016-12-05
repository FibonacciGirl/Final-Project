

n=120
data<-array(dim=c(n,n.dimensions,n.))

ny=60
nx=60
n.dimensions=2

y<-matrix(nrow=ny,ncol=n.dimensions)
x<-matrix(nrow=nx,ncol=n.dimensions)

y[(1:30),1]<-rep(seq(7,16,by=1),3)
y[(1:30),2]<-c(rep(1,10),rep(2,10),rep(3,10))
y[(31:60),1]<-c(rep(1,10),rep(2,10),rep(3,10))
y[(31:60),2]<-rep(seq(from=7,to=16,by=1),3)


x[(1:30),1]<-rep(seq(7,16,by=1),3)
x[(1:30),2]<-c(rep(4,10),rep(5,10),rep(6,10))
x[(31:60),1]<-c(rep(4,10),rep(5,10),rep(6,10))
x[(31:60),2]<-rep(seq(from=7,to=16,by=1),3)



print(x)
mfrow(c(1,2))
plot(y)
plot(x)
for(s in 1:n){
  for(i to 1:n.dimensions){
    X[s,i]<-
  }
}

x[,1]

x[1:30,1]<-seq(1:10)
x[31:60,1]<-seq(4:6)
x[1:30,2]<-seq(1:3)
x[31:60,2]<-seq(1:10)

X.dat<-sample(x,60)
Y.dat<-sample(y,60)
cat<-c(rep(1,60),rep(2,60))


data<-rbind(X.dat,Y.dat,cat)
data
plot(x,y)
plot(y)


training.data<-data.matrix(data[,sample(1:120,120)])
training.data


install.packages('ggplot2')
library(ggplot2)
ggplot(training.data, mapping=aes(x=training.data[1,],y=training.data[2,],category=))
