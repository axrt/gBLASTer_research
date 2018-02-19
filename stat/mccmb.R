#load the data
master.table.5<-read.table(file = "master_table5_1.txt",header = TRUE,sep = "\t",row.names = 1)

#select 3,6
master.table.5.3x6<-master.table.5[grepl(pattern = "\\b[3,6]X",x = rownames(master.table.5),perl = TRUE),c(-3,-6)]
master.table.5.3x6.hclust<-hclust(d = dist(master.table.5.3x6,method = "eu"),method = "ward.D")
plot(master.table.5.3x6.hclust)

#cluster
library("gbra")
master.table.5.3x6.hclust.cluster<-data.frame(order=rownames(master.table.5.3x6)[master.table.5.3x6.hclust$order])
write.table(x = master.table.5.3x6.hclust.cluster,file = "master_table5_1_3X6_nocol.cluster",sep = "\t",quote = FALSE)
dnd<-color.clust(gen.hclust = master.table.5.3x6.hclust,colors = c("grey","black"),sep = "X",g.ids = c(3,6))

pdf(file = "master_table5_1_3X6_nocol.pdf",width = 25,height = 25)
plot(dnd,leaflab = "none")
dev.off()

#logistic regression
attach(master.table.5.3x6)
master.table.5.3x6$cluster<-grepl(pattern = "\\b[6]X",x = rownames(master.table.5.3x6),perl = TRUE)
glm.fit<-glm(formula = as.formula(paste("cluster ~",paste(colnames(master.table.5.3x6[1:ncol(master.table.5.3x6)-1]),collapse = "+"))),data=master.table.5.3x6,family=binomial)
summary(glm.fit)
glm.fit<-glm(formula = as.formula("cluster ~ X1+X4+X5+X10+X15+X17+X22+X23+X34+X37+X46"),data=master.table.5.3x6,family=binomial)
summary(glm.fit)
glm.probs<-predict(glm.fit,type = "response")
glm.probs.drosoph<-sapply(glm.probs,function(i){if(i>=0.5){return(TRUE)}else{return(FALSE)}})
table(glm.probs.drosoph,cluster)
cross<-data.frame(cluster=cluster,probs=glm.probs)
cross$match<-apply(cross,1,function(i){
  return(i[1]==(i[2]>0.5))
})
cross.mismatch<-cross[cross$match==FALSE,]
cross.mismatch[cross.mismatch$probs<=0.05,]
nrow(cross.mismatch[cross.mismatch$probs>=0.95,])
summary(cross$probs)
write.table(x = cross.mismatch,file = "cross.mismatch.master.5.3X6.txt",quote = FALSE,sep = "\t",col.names = colnames(cross.mismatch))


#QDA
require("MASS")
master.table.5.3x6<-master.table.5[grepl(pattern = "\\b[3,6]X",x = rownames(master.table.5),perl = TRUE),c(-3,-6)]
master.table.5.3x6$cluster<-grepl(pattern = "\\b[6]X",x = rownames(master.table.5.3x6),perl = TRUE)
#master.table.5.3x6<-as.data.frame(apply(master.table.5.3x6,2,jitter))
master.table.5.3x6$cluster<-sapply(rownames(master.table.5.3x6),function(i){if(grepl(pattern = "\\b[6]X",x = i,perl = TRUE)){
  return(TRUE)
}else{
  return(FALSE)
}})
qda.fit<-qda(formula = as.formula(paste("cluster ~ X1+X4+X5+X10+X15+X17+X22+X23+X34+X37+X46")),data=master.table.5.3x6,subet=train)
qda.probs<-predict(qda.fit,master.table.5.3x6)
qda.qal<-table(qda.probs$class,cluster)
qda.true.positive<-(qda.qal[1,1]+qda.qal[2,2])/sum(qda.qal)
qda.cross<-data.frame(cluster=cluster,probs=qda.probs$posterior)
head(qda.cross)
qda.cross$match<-sapply(1:nrow(qda.cross),function(i){
  if(grepl(pattern = "\\b[3]X",x =row.names(qda.cross)[i],perl = TRUE)){
    return(qda.cross[i,2]>=0.5)
  }else{
    return(qda.cross[i,2]<0.5)
  }
})
qda.cross.mismatch<-qda.cross[qda.cross$match==FALSE,]
length(intersect(rownames(cross.mismatch),rownames(qda.cross.mismatch)))

#LDA
master.table.5.3x6<-master.table.5[grepl(pattern = "\\b[3,6]X",x = rownames(master.table.5),perl = TRUE),c(-3,-6)]
master.table.5.3x6$cluster<-grepl(pattern = "\\b[6]X",x = rownames(master.table.5.3x6),perl = TRUE)
lda.fit<-lda(formula = as.formula(paste("cluster ~ X1+X4+X5+X10+X15+X17+X22+X23+X34+X37+X46")),data=master.table.5.3x6,subet=train)
lda.probs<-predict(lda.fit,master.table.5.3x6)
lda.pred<-sapply(lda.probs$posterior[1:nrow(master.table.5.3x6)],function(i){if(i>=0.5){return(FALSE)}else{return(TRUE)}})
lda.qual<-table(lda.pred,cluster)
lda.qual<-(lda.qual[1,1]+lda.qual[2,2])/sum(lda.qual)

lda.missmatch<-rownames(master.table.5.3x6)[lda.pred!=cluster]

length(intersect(rownames(cross.mismatch),intersect(rownames(qda.cross.mismatch),lda.missmatch)))

