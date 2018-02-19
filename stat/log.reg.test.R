library(gbra)
library(doParallel)
library(reshape2)
library(ggbiplot)
bh.data<-read.bhs(bh.folder = "/home/alext/Documents/Research/gBLASTer/bh",sep="_") %>% 
  as.data.frame() %>%
  restrict.minimal.hits(minhit = 4) %>%
  normalize.scores() %>%
  attach.genomeid.header() %>%
  sign.bh.table()
write.table(x=bh.data, file = "bh.data.minhit_4.txt", sep="\t", quote = FALSE)
legend<-load.legend(legend.file = "/home/alext/Developer/gbra/inst/extradata/legend.csv",sep = "\t")

############# 3 vs 14 ############# 
org1=3
org2=14
data.pair<-bh.data %>% select.genomes(g.ids = c(org1,org2))
data.pair$cluster<-grepl(pattern = paste("\\b",org2,"X",sep = "",collapse = ""),x = rownames(data.pair),perl = TRUE)
table(data.pair$cluster)
model.3.14<-logreg.mismatch.genes(data = data.pair,org1 = org1,org2 =org2)
summary(model.3.14$fit)

data.pair.stn<-as.data.frame(dist.stdev.one(df = data.pair))
model.3.14.stn<-logreg.mismatch.genes(data = data.pair.stn,org1 = org1,org2 =org2)
summary(model.3.14.stn)

model.3.14$mismatch
model.3.14.stn$mismatch
intersect(rownames(model.3.14$mismatch),rownames(model.3.14.stn$mismatch))

############# 3 vs 4 ############# 
org1=3
org2=4
data.pair<-bh.data %>% select.genomes(g.ids = c(org1,org2))
data.pair$cluster<-grepl(pattern = paste("\\b",org2,"X",sep = "",collapse = ""),x = rownames(data.pair),perl = TRUE)
table(data.pair$cluster)
model.3.4<-logreg.mismatch.genes(data = data.pair,org1 = org1,org2 =org2)
summary(model.3.4$fit)

data.pair.stn<-as.data.frame(dist.stdev.one(df = data.pair))
model.3.4.stn<-logreg.mismatch.genes(data = data.pair.stn,org1 = org1,org2 =org2)
summary(model.3.4.stn$fit)

model.3.4$mismatch
model.3.4.stn$mismatch
intersect(rownames(model.3.4$mismatch),rownames(model.3.4.stn$mismatch))

############# 3 vs 6 ############# 
org1=3
org2=6
data.pair<-bh.data %>% select.genomes(g.ids = c(org1,org2))
data.pair$cluster<-grepl(pattern = paste("\\b",org2,"X",sep = "",collapse = ""),x = rownames(data.pair),perl = TRUE)
table(data.pair$cluster)
model.3.6<-logreg.mismatch.genes(data = data.pair,org1 = org1,org2 =org2)
summary(model.3.6$fit)

data.pair.stn<-as.data.frame(dist.stdev.one(df = data.pair))
model.3.6.stn<-logreg.mismatch.genes(data = data.pair.stn,org1 = org1,org2 =org2)
summary(model.3.6.stn$fit)

model.3.6$mismatch
model.3.6.stn$mismatch
intersect(rownames(model.3.6$mismatch),rownames(model.3.6.stn$mismatch))


############# 1 vs 2 ############# 
org1=1
org2=2
data.pair<-bh.data %>% select.genomes(g.ids = c(org1,org2))
data.pair$cluster<-grepl(pattern = paste("\\b",org2,"X",sep = "",collapse = ""),x = rownames(data.pair),perl = TRUE)
table(data.pair$cluster)
model.1.2<-logreg.mismatch.genes(data = data.pair,org1 = org1,org2 =org2)
summary(model.1.2$fit)

data.pair.stn<-as.data.frame(dist.stdev.one(df = data.pair))
model.1.2.stn<-logreg.mismatch.genes(data = data.pair.stn,org1 = org1,org2 =org2)
summary(model.1.2.stn$fit)

model.1.2$mismatch
model.1.2.stn$mismatch
intersect(rownames(model.1.2$mismatch),rownames(model.1.2.stn$mismatch))


############# 14 vs 44 ############# 
org1=44
org2=14
data.pair<-bh.data %>% select.genomes(g.ids = c(org1,org2))
data.pair$cluster<-grepl(pattern = paste("\\b",org2,"X",sep = "",collapse = ""),x = rownames(data.pair),perl = TRUE)
table(data.pair$cluster)
model.14.44<-logreg.mismatch.genes(data = data.pair,org1 = org1,org2 =org2)
summary(model.14.44$fit)

data.pair.stn<-as.data.frame(dist.stdev.one(df = data.pair))
model.14.44.stn<-logreg.mismatch.genes(data = data.pair.stn,org1 = org1,org2 =org2)
summary(model.14.44.stn$fit)

model.14.44$mismatch
model.14.44.stn$mismatch
intersect(rownames(model.14.44$mismatch),rownames(model.14.44.stn$mismatch))


############# 7 vs 17 ############# 
org1=7
org2=17
data.pair<-bh.data %>% select.genomes(g.ids = c(org1,org2))
data.pair$cluster<-grepl(pattern = paste("\\b",org2,"X",sep = "",collapse = ""),x = rownames(data.pair),perl = TRUE)
table(data.pair$cluster)
model.7.17<-logreg.mismatch.genes(data = data.pair,org1 = org1,org2 =org2)
summary(model.7.17$fit)

data.pair.stn<-as.data.frame(dist.stdev.one(df = data.pair))
model.7.17.stn<-logreg.mismatch.genes(data = data.pair.stn,org1 = org1,org2 =org2)
summary(model.7.17.stn$fit)

model.7.17$mismatch
model.7.17.stn$mismatch
intersect(rownames(model.7.17$mismatch),rownames(model.7.17.stn$mismatch))


############# 45 vs 16 ############# 
org1=45
org2=16
data.pair<-bh.data %>% select.genomes(g.ids = c(org1,org2))
data.pair$cluster<-grepl(pattern = paste("\\b",org2,"X",sep = "",collapse = ""),x = rownames(data.pair),perl = TRUE)
table(data.pair$cluster)
model.45.16<-logreg.mismatch.genes(data = data.pair,org1 = org1,org2 =org2)
summary(model.45.16$fit)

data.pair.stn<-as.data.frame(dist.stdev.one(df = data.pair))
model.45.16.stn<-logreg.mismatch.genes(data = data.pair.stn,org1 = org1,org2 =org2)
summary(model.45.16.stn$fit)

model.45.16$mismatch
model.45.16.stn$mismatch
intersect(rownames(model.45.16$mismatch),rownames(model.45.16.stn$mismatch))

############# 9 vs 17 ############# 
org1=9
org2=17
data.pair<-bh.data %>% select.genomes(g.ids = c(org1,org2))
data.pair$cluster<-grepl(pattern = paste("\\b",org2,"X",sep = "",collapse = ""),x = rownames(data.pair),perl = TRUE)
table(data.pair$cluster)
model.9.17<-logreg.mismatch.genes(data = data.pair,org1 = org1,org2 =org2)
summary(model.9.17$fit)

data.pair.stn<-as.data.frame(dist.stdev.one(df = data.pair))
model.9.17.stn<-logreg.mismatch.genes(data = data.pair.stn,org1 = org1,org2 =org2)
summary(model.9.17.stn$fit)

model.9.17$mismatch
model.9.17.stn$mismatch
intersect(rownames(model.45.16$mismatch),rownames(model.9.17.stn$mismatch))

###################################

total.model.list<-total.logit.model(df = bh.data, legend = legend[1:(nrow(legend)),],cores = 12)

###################################

data.pca<-data.pca %>% select.genomes(df = ., g.ids = c(3,14))
data.pca.stn<-dist.stdev.one(df=data.pca)
gen.pca(df = data.pca, legend = legend, g.ids = c(3,14))
gen.pca(df = data.pca.stn, legend = legend, g.ids = c(3,14))

###################################

data.rsq<-data.frame(matrix(ncol=nrow(legend),nrow=nrow(legend),1),row.names = legend$name)
colnames(data.rsq)<-legend$name
dir<-"data"
for(i in 1:(nrow(legend)-1)){
  for(j in (i+1):nrow(legend)){
    if(i!=j){
      model<-load(file = paste(dir,paste(paste(legend$id_genomes[i],legend$id_genomes[j],sep="_"),".rda",sep = ""),sep = "/"))
      rs<-as.numeric(l.m.g$Rsq)
      data.rsq[i,j]<-rs
      data.rsq[j,i]<-rs
      message(rs)
    }
  }
}
write.table(x=data.rsq, file="data.rsq.txt", sep="\t")

data.tp<-data.frame(matrix(ncol=nrow(legend),nrow=nrow(legend),0),row.names = legend$name)
colnames(data.tp)<-legend$name
counter<-0
for(i in 1:(nrow(legend)-1)){
  for(j in (i+1):nrow(legend)){
    if(i!=j){
      model<-load(file = paste(dir,paste(paste(legend$id_genomes[i],legend$id_genomes[j],sep="_"),".rda",sep = ""),sep = "/"))
      tab<-model.ratio.table(model = l.m.g$fit, predict.data = l.m.g$working.data)
      tp<-0
      if(nrow(tab)>1){
        tp<-(tab[1,1]+tab[2,2])/sum(tab)
      }
      data.tp[i,j]<-tp
      data.tp[j,i]<-tp
      counter<<-counter+1
      message(paste(counter/((nrow(legend)^2)/2-nrow(legend))*100,"% done with ",tp," value",sep = ""))
    }
  }
}
write.table(x=data.tp, file="data.tp.txt", sep="\t")


data.pval<-data.frame(matrix(ncol=nrow(legend),nrow=nrow(legend),0),row.names = legend$name)
colnames(data.pval)<-legend$name
counter<-0
for(i in 1:(nrow(legend)-1)){
  for(j in (i+1):nrow(legend)){
    if(i!=j){
      model<-load(file = paste(dir,paste(paste(legend$id_genomes[i],legend$id_genomes[j],sep="_"),".rda",sep = ""),sep = "/"))
      s<-summary(l.m.g$fit)
      min.pv<-min(s$coefficients[,4])
      data.pval[i,j]<-min.pv
      data.pval[j,i]<-min.pv
      counter<<-counter+1
      message(paste(counter/((nrow(legend)^2)/2-nrow(legend))*100,"% done with ",min.pv," value",sep = ""))
    }
  }
}
write.table(x=data.pval, file="data.pval.txt", sep="\t")

library("reshape2")
data.rsq<-read.table(file = "/home/alext/Documents/Research/gBLASTer/stat/data.rsq.txt", stringsAsFactors = FALSE)
for(i in 1:(nrow(legend))){
  for(j in (i):nrow(legend)){
    if(i==j){
      data.rsq[i,j]<-0
    }
  }
}
data.tp.melt<-melt(data.tp)
data.rsq.melt<-melt(data.rsq)
data.melt<-cbind(data.tp.melt, data.rsq.melt)[,c(1,2,4)]
colnames(data.melt)<-c("organism", "tp", "rsq")

gp<-ggplot(data = data.melt, mapping=aes(x=log(tp), y=log(rsq), color=organism))+geom_line()+facet_grid(organism~.)+facet_wrap(~organism,ncol = 5)
pdf(file="rsq.tp.pdf", width=25, height = 25)
plot(gp)
dev.off()
data.melt$ltp<-log(data.melt$tp)
data.melt$lrsq<-log(data.melt$rsq)
data.melt.pure<-data.melt %>% filter(is.finite(ltp), is.finite(lrsq))


tp.rsq.model<-lm(data.melt.pure,formula=ltp~lrsq)
summary(tp.rsq.model)
################################### for dr. Syvanen
data.3.6 <- select.genomes(df = bh.data, g.ids = c(3,6))
pc<-gen.pca(df = data.3.6, g.ids = c(3,6),legend = legend)
png(filename = "plot1.png",width=700, height = 700)
plot(pc)
dev.off()

data.3.6<-data.3.6[,!colnames(data.3.6)%in%c("X3","X6")]

pc2<-color.clust(gen.hclust = hclust(dist(data.3.6),method = "ward.D2"),g.ids = c(3,6),colors = c("red","blue"))
png(filename = "plot2.png",width=700, height = 700)
plot(pc2,leaflab="none")
dev.off()

data.2.38 <- select.genomes(df = bh.data, g.ids = c(2,38))
data.2.38<-data.2.38[,!colnames(data.2.38)%in%c("X2","X38")]
pc<-gen.pca(df = data.2.38, g.ids = c(2,38),legend = legend)
plot(pc)
pc2<-color.clust(gen.hclust = hclust(dist(data.2.38),method = "ward.D2"),g.ids = c(2,38),colors = c("red","blue"))
plot(pc2,leaflab="none")

data.1.3 <- select.genomes(df = bh.data, g.ids = c(1,3))
data.1.3<-data.1.3[,!colnames(data.1.3)%in%c("X1","X3")]
pc<-gen.pca(df = data.1.3, g.ids = c(1,3),legend = legend)
plot(pc)
pc2<-color.clust(gen.hclust = hclust(dist(data.1.3),method = "ward.D2"),g.ids = c(1,3),colors = c("red","blue"))
plot(pc2,leaflab="none")


data.1.2 <- select.genomes(df = bh.data, g.ids = c(1,2))
data.1.2<-data.1.2[,!colnames(data.1.2)%in%c("X1","X2")]
pc<-gen.pca(df = data.1.2, g.ids = c(1,2),legend = legend)
plot(pc)
pc2<-color.clust(gen.hclust = hclust(dist(data.1.2),method = "ward.D2"),g.ids = c(1,2),colors = c("red","blue"))
plot(pc2,leaflab="none")
###################################

factor.core<-c("core","outsider","nosep")
master.table.core<-as.data.frame(matrix(nrow=nrow(bh.data),ncol=nrow(legend),factor.core[3]),stringsAsFactors = FALSE)
row.names(master.table.core)<-row.names(bh.data)
colnames(master.table.core)<-legend$id_genomes

for(i in 1:(nrow(legend)-1)){
  for(j in (i+1):nrow(legend)){
    print(paste(i,j))
    if(i!=j){
      load(file = paste(dir,paste(paste(legend$id_genomes[i],legend$id_genomes[j],sep="_"),".rda",sep = ""),sep = "/"))
      if(as.numeric(l.m.g$Rsq)>=0.75){
        data.table<-select.genomes(df = bh.data, g.ids = c(legend$id_genomes[i],legend$id_genomes[j]))
        glm.probs<-predict(l.m.g$fit, data.table, type = "response")
        glm.probs.groups<-sapply(glm.probs,function(i){if(i>=0.5){return(FALSE)}else{return(TRUE)}})
        data.table$cluster<-grepl(pattern = paste("\\b",legend$id_genomes[i],"X",sep = "",collapse = ""),x = rownames(data.table),perl = TRUE)
        sep.factor.v<-ifelse(test = data.table$cluster==glm.probs.groups, yes = factor.core[1], no = factor.core[2])
        sep.factor.v.i<-sep.factor.v[grepl(pattern = paste("\\b",legend$id_genomes[i],"X",sep = "",collapse = ""),x = names(sep.factor.v),perl = TRUE)]
        sep.factor.v.j<-sep.factor.v[grepl(pattern = paste("\\b",legend$id_genomes[j],"X",sep = "",collapse = ""),x = names(sep.factor.v),perl = TRUE)]
        
        core.col.j<-which(legend$id_genomes==legend$id_genomes[j])
        core.col.i<-which(legend$id_genomes==legend$id_genomes[i])
        for(k in 1:length(sep.factor.v.i)){
          master.table.core[names(sep.factor.v.i)[k],core.col.j]<-sep.factor.v.i[k]
        }
        for(k in 1:length(sep.factor.v.j)){
          master.table.core[names(sep.factor.v.j)[k],core.col.i]<-sep.factor.v.j[k]
        }
      }
    }
  }
}
write.table(x=master.table.core, file="master.table.core.txt", sep="\t", quote = FALSE)

master.table.core<-read.table(file="master.table.core.txt")
cluster <<- makeCluster(6,outfile = "cluster.txt")
registerDoParallel(cluster)
master.table.core.factor<-parApply(cl = cluster, X = master.table.core, MARGIN = 1, FUN = function(i){
  return(sapply(i, function(x){
    if(x=="core"){return(0)}
    if(x=="outsider"){return(1)}
    if(x=="nosep"){return(2)}
  }))
})
stopCluster(cluster)
master.table.core.factor<-t(master.table.core.factor)
pca.core<-prcomp(master.table.core.factor)
plot(pca.core)
pca.core.data<-cbind(pca.core$x[,1],pca.core$x[,2])
colnames(pca.core.data)<-c("PC1", "PC2")
pca.core.data<-as.data.frame(pca.core.data)
groups<-sapply(1:nrow(master.table.core.factor), function(i){return(strsplit(rownames(master.table.core.factor)[i],"X")[[1]][1])})
pca.core.data$groups<-legend$name
core.plot<-ggplot(data=pca.core.data, mapping=aes(x=PC1, y=PC2, color=groups))+geom_point(aes(alpha=0.1))
core.plot
pdf(file = "core.plot.pdf", width=10, height=10)
ggbiplot(pca.core, obs.scale = 1, var.scale = 1, groups = groups, alpha = 0.1)
dev.off()

master.table.core.summary<-t(master.table.core) %>% melt(.,id.vars = colnames(.)) %>% group_by(Var1,value) %>%  
  summarise(sum=n()) %>% group_by(Var1) %>% do({
    return(cbind(.,ratio=rescale(.$sum,from=c(0,sum(.$sum)), to=c(0,1))))
  })

master.table.core.summary.core.outsider<- master.table.core.summary %>% filter(value!="nosep")
cluster <<- makeCluster(6,outfile = "cluster.txt")
clusterExport(cluster,"master.table.core.summary.core.outsider")
registerDoParallel(cluster)
master.table.core.summary.core.outsider$genome<-parSapply(cl = cluster, X = 1:nrow(master.table.core.summary.core.outsider), FUN = function(i){return(strsplit(rownames(master.table.core.summary.core.outsider)[i],"X")[[1]][1])})
stopCluster(cluster)
pdf(file="core.ratio.pdf",width = 25, height=5)
ggplot(data=master.table.core.summary, mapping=aes(x=value, y=ratio, fill=Var1))+geom_bar(stat="identity")+facet_grid(.~Var1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),legend.position="top")
dev.off()

master.table.core.plot<-t(master.table.core) %>% melt(.,id.vars = colnames(.)) %>% group_by(Var2,value) %>%  
  summarise(sum=n())
cluster <<- makeCluster(6,outfile = "cluster.txt")
clusterExport(cluster,"master.table.core.plot")
registerDoParallel(cluster)
master.table.core.plot$genome<-parSapply(cl = cluster, X = 1:nrow(master.table.core.plot), 
                   FUN = function(i){return(strsplit(as.character(master.table.core.plot$Var2[i]),"X")[[1]][1])})
stopCluster(cluster)

curr<-0
master.table.core.plot<-master.table.core.plot %>% group_by(Var2) %>% do({
  if(.[1,4]!=curr){
    curr<<-.[1,4]
    print(paste("Genome changed to:",.[1,4]))
  }
  if(as.character(.[nrow(.),1])=="46X17701048"){print("done")}
  return(data.frame(., ratio=sapply(.[,3],function(i){return(i/45*100)} )))
}) %>% group_by(value) %>% arrange(ratio)


ggplot(data=., mapping=aes(x=Var2,y=sum, color=ratio))+geom_point()+facet_grid(.~Var1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),legend.position="top")

master.table.core.plot.46<-master.table.core.plot %>% filter(genome==46)
master.table.core.plot.4<-master.table.core.plot %>% filter(genome==4)%>% do({
  if(.[1,4]!=curr){
    curr<<-.[1,4]
    print(paste("Genome changed to:",.[1,4]))
  }
  if(as.character(.[nrow(.),1])=="46X17701048"){print("done")}
  return(data.frame(., ratio=sapply(.[,3],function(i){return(i/45*100)} )))
})


master.table.core.plot$genome<-factor(master.table.core.plot$genome)
master.table.core.plot<-master.table.core.plot %>% group_by(Var2) %>% group_by(value) %>% arrange(sum)
master.table.core.plot$Var2<-factor(master.table.core.plot$Var2, levels=unique(master.table.core.plot$Var2))
master.table.core.plot %>% group_by(genome) %>% do({
  pdf(file=paste("each.gene.",.[1,4],".pdf",sep=""), width=0.05*nrow(.), height=6)
  plot(ggplot(data=., mapping=aes(x=Var2, y=sum, group=value, fill=value, color=value))+geom_line(alpha=0.333)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),legend.position="top"))
  dev.off()
  return(data.frame())
})

master.table.core.plot$ratio<-master.table.core.plot$sum/45
master.table.core.plot$genome<-as.numeric(as.character(master.table.core.plot$genome))
  
gp<-master.table.core.plot %>% arrange(genome) %>% ggplot(data=., mapping=aes(x=ratio, fill=value,color=value))+
  geom_histogram(binwidth=0.2, position = "dodge")+facet_grid(~genome)+ facet_wrap(~genome,ncol=9)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),legend.position="top")

pdf(file="each.gene.pdf", width=25, height=10)
plot(gp)
dev.off()

gp<-master.table.core.plot %>% arrange(genome) %>% ggplot(data=.)+
  geom_histogram(aes(x=ratio, y=..ncount.., color=value, fill=value),binwidth=0.01, position = "dodge")+facet_grid(~genome)+ facet_wrap(~genome,ncol=9)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),legend.position="top")

pdf(file="each.gene.normal.pdf", width=25, height=10)
plot(gp)
dev.off()


master.table.core.cut_3<-master.table.core.plot %>% filter(value=="core", sum>=3)
write.table(x=master.table.core.cut_3, file="master.table.core.cut_3.txt",sep = "\t", quote = FALSE)
master.table.outsider.cut_3<-master.table.core.plot %>% filter(value=="outsider", sum>=3)
write.table(x=master.table.outsider.cut_3, file="master.table.outsider.cut_3.txt",sep = "\t", quote = FALSE)
