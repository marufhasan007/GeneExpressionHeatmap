library(openxlsx)
setwd("Z:/30Deser/user/hasan/Privat/PCR_Phosphorus_Low_High/R")
outlier<-"n"     # y or n
dat<-read.csv("Z:/30Deser/user/hasan/Privat/PCR_Phosphorus_Low_High/LightCycler480/qPCR_Result.csv") #edited file
dat<-dat[,c(1,5,6,7)]
dat$gene<-sub("_.*","",dat$Experiment.Name)

dat<-dat[!dat$SampleName=="NC",]
dat<-dat[-grep("std.*",dat$SampleName),]
dat<-dat[,c(5,2,3,4)]
colnames(dat)<-c("gene","Sample","CP","Conc")
data<-dat[,-3]
data<-data[order(data[,1],data[,2]),]

data1 <- data[seq_along(data[,1]) %% 2==0,]
data1$ID<-paste0(data1$gene,"x",data1$Sample)
data2 <- data[!seq_along(data[,1]) %% 2==0,]
data2$ID<-paste0(data2$gene,"x",data2$Sample)
dat.comb<-merge(data1,data2,by="ID")
dat.comb<-dat.comb[,c(3,2,4,7)]
colnames(dat.comb)<-c("Sample","gene","rep1","rep2")

a<-split(dat.comb,dat.comb$gene)
a<-lapply(a, function (x) {colnames(x)[c(3,4)]<-c(paste0(x[1,2],1),paste0(x[1,2],2));x})
a<-lapply(a, function(x) { x["gene"] <- NULL; x })
res<-Reduce(function(x, y) merge(x, y, all=TRUE,by="Sample"), a)

pheno<-cbind.data.frame(Sample=res$Sample,tissue=gsub("[0-9]+","",res$Sample),ID=sub("[A-z]+","",res$Sample))
dat<-merge(pheno,res,by.x="Sample",by.y="Sample")

##############################
###          Edit!!!       ###
##############################
datum <- c("2708021");  #date
housekeeping_genes <- c("RPL32")
ErstesTranscript <- "CYP24A1";   ### give name of first gene (excluding HK)
norm_method <- "per_group"       ### "per_group", or "all"
normgroup   <- "tissue"          ### if "per_group" provide the group
###############################

### reorder columns that HK come first ###
Tieranzahl = dim(dat)[1]; #number of samples
pos_HK<-grep(housekeeping_genes[1],colnames(dat))
pos_gene<-grep(ErstesTranscript[1],colnames(dat))
if (pos_HK[1]<pos_gene[1]){
  dat<-dat[,c(1:(min(pos_HK)-1),pos_HK,(max(pos_HK)+1):dim(dat)[2])]    
} else
  dat<-dat[,c(1:(min(pos_gene)-1),pos_HK,pos_gene,(max(pos_gene)+1):(min(pos_HK)-1),(max(pos_HK)+1):dim(dat)[2])]
colnames(dat)<-gsub("[0-9]$","",colnames(dat))

housekeeping_genes_out <- c()
if(length(housekeeping_genes_out)>0){
  dat <- dat[,-grep(housekeeping_genes_out, colnames(dat))]
} else
  rm(housekeeping_genes_out)

pos <- match(housekeeping_genes[1], colnames(dat))                              

#outlier test
library(outliers)
library(ggplot2)
grubbs.flag <- function(x) {
  outliers <- NULL
  test <- x
  grubbs.result <- grubbs.test(test)
  pv <- grubbs.result$p.value
  while(pv < 0.05) {
    outliers <- c(outliers,as.numeric(strsplit(grubbs.result$alternative," ")[[1]][3]))
    test <- x[!x %in% outliers]
    grubbs.result <- grubbs.test(test)
    pv <- grubbs.result$p.value
  }
  return(data.frame(X=x,Outlier=(x %in% outliers)))
}

dat.outtest <- dat[,c(pos:dim(dat)[2])]
dat.outtest2a <- dat.outtest[,seq(1,by=2,ncol(dat.outtest))]
dat.outtest2b <- dat.outtest[,seq(2,by=2,ncol(dat.outtest))]
colnames(dat.outtest2b) <- colnames(dat.outtest2a)
dat.outtest <- rbind(dat.outtest2a, dat.outtest2b)

dat.outtest.list <- apply(dat.outtest, 2, grubbs.flag)

dat.out2 <- data.frame(ID=c(1:Tieranzahl))
for(i in c(1:dim(dat.outtest)[2])) {
  truefalse <- as.data.frame(dat.outtest.list[i])[,2]
  dat.out2 <- cbind(dat.out2, truefalse)
  colnames(dat.out2)[i+1] <- colnames(dat.outtest)[i]
}
dat.out2 <- cbind(dat.out2[1:(nrow(dat.out2)/2),],dat.out2[-(1:(nrow(dat.out2)/2)),])[,c(matrix(data = 1:(ncol(dat.out2)*2), nrow=2, byrow=TRUE))]
dat.out2 <- dat.out2[,-2]
#dat.out2

#write.table(dat.out2, paste(datum,"qPCR_outlier.txt",sep="_"), col.names=TRUE, row.names=FALSE, quote=FALSE,sep='\t')

#### remove outlier #####
# no outlier found for u, we skip this part ... 
if (outlier=="y"){ # something is missing here ... outlier is not defined correctly
  df2<-cbind("ID"=dat.out2[,1],dat[,c(pos:dim(dat)[2])])
  remov<-which(dat.out2=="TRUE",arr.ind=T)
  
  if (dim(remov)[1]>0) {
    for (i in 1:dim(remov)[1])
      df2[remov[i,1],remov[i,2]]<-NA
    dat<-cbind(dat[,c(1:pos-1)],df2[,c(2:dim(df2)[2])])
  }else { dat<-dat}
}

######

#Generate means per transcript from replicates
dat2 <- data.frame(ID=c(1:Tieranzahl))

for(i in c(1: c(length(pos:dim(dat)[2])/2) )){
  Means <- rowMeans(dat[,(2*i+pos-2):(2*i+pos-1)],na.rm=TRUE)
  dat2  <- cbind(dat2, Means)
  colnames(dat2)[i+1] <- colnames(dat[(2*i+pos-2)])
} 
dat2  <- cbind(dat[,(1:c(pos-1))], dat2[,2:dim(dat2)[2]])
ui<-sapply(dat2,is.nan)
dat2[ui]<-NA
#dat2


#mean per group or over all samples using housekeeping genes
dat3 <- dat2
if(norm_method == "per_group") {              
  for(i in match(housekeeping_genes, colnames(dat2))){                    # Mittelwertberechnungen per Gruppe f?r housekeeping_genes
    means=sapply(split(dat3[1:Tieranzahl,i], dat3[1:Tieranzahl, match(normgroup, colnames(dat2))]), mean, na.rm=TRUE)
    dat3[,i]=means[dat3[,match(normgroup, colnames(dat2))]]
  }
} else
  for(i in match(housekeeping_genes, colnames(dat2))){                    # Mittelwertberechnungen ?ber alle f?r housekeeping_genes
    means=mean(dat3[1:Tieranzahl,i],na.rm=TRUE)
    dat3[,i]=means
  }
#dat3

#geometric mean of averaged values based on housekeeping genes
norm_fac1 <- dat2[,match(housekeeping_genes, colnames(dat2))]/dat3[,match(housekeeping_genes, colnames(dat3))]                        
norm_fac2 <- c(1:Tieranzahl)
if(length(housekeeping_genes)>1){
  norm_fac2=apply(norm_fac1, 1, function(x) exp(mean(log(x),na.rm=TRUE)))             #geometric.mean(norm_fac[i,]) per row
}else
  norm_fac2=norm_fac1        

#divid target genes by the geometric mean of  housekeeping_genes, rename +combine +transpose colnames/rownames
qPCR <- dat2[,match(ErstesTranscript, colnames(dat2)):dim(dat2)[2]] / norm_fac2
colnames(qPCR) <- paste(colnames(qPCR),"qPCR",sep = "_")
rownames(qPCR) <- dat2[,1]
q.stat<-qPCR
qPCR <- as.data.frame(t(qPCR))
qPCR <- cbind(sub("_qPCR", "", rownames(qPCR)), qPCR)
colnames(qPCR)[1] <- "gene"
write.csv(qPCR, paste0(paste(datum,"qPCR_VitD", norm_method, normgroup, paste(housekeeping_genes, collapse="_"),sep="_"),".csv"))
qPCR
write.xlsx (qPCR, "Data_Normalized.xlsx" )