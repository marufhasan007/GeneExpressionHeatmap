setwd("Z:/30Deser/user/hasan/Privat/PCR_Phosphorus_Low_High/R")
x <- read.csv("Z:/30Deser/user/hasan/Privat/PCR_Phosphorus_Low_High/R/PCR_Low.csv")
rownames(x)<-x[,1] # First row is labeled as the name of the row
x<-x[,-1] # To delete the repeated column you will see on console

### mean, with NA or 0?

tis<-unique(sub("[0-9]+","",colnames(x)))  
a <- matrix(nrow = dim(x)[1], ncol = length(tis))
colnames(a)<-tis
rownames(a)<-rownames(x)

for (i in tis){
  a[,i]<-rowMeans(x[,grep(i,colnames(x))],na.rm=TRUE)
}
alog<-log(a,2)
write.csv(alog,"PCR_Low.csv")


