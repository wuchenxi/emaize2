library(rrBLUP)
library(data.table)

ntrain <- 4754

T <- fread('chr1_emaize.genoMat',nrows=20)
P <- fread('pheno_emaize.txt')
nr <- nrow(T)
nc <- ncol(T)

X<-T[,5:nc]
Y<-data.matrix(P[1:ntrain,4])
for(i in 1:nr){
   symb<-T[[i,2]]
   
   major<-substring(symb,1,1)
   minor<-substring(symb,3,3)
   val<-c(-1,0,0,1)
   names(val)<-c(paste(major,major,sep=""),paste(major,minor,sep=""),paste(minor,major,sep=""),paste(minor,minor,sep=""))
   X[i]<-lapply(X[i],function(x){val[[x]]})
}
Xtrain<-data.matrix(X[,1:ntrain])
#colnames(Xtrain)<-1:ntrain
#geno<-data.frame(marker=T[,1],chrom=T[,3],pos=T[,4],Xtrain,check.names=FALSE)
#pheno<-data.frame(line=1:ntrain,y=Y)

mixed.solve(y=data.matrix(P[,4]),K=A.mat(t(data.matrix(X))),X=t(data.matrix(X)))