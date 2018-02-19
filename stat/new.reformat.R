library(gbra)
master.table.raw<-read.bhs(bh.folder = "/home/alext/Documents/Research/gBLASTer/bh")
lcr.acs<-as.numeric(read.table(file = "/home/alext/Documents/Research/gBLASTer/orfs/lcrs/orfs.lcracs",header = FALSE)$V1)
length(setdiff(master.table.raw$QUERY_ORF_ID, lcr.acs))
master.table.raw<-master.table.raw[master.table.raw$QUERY_ORF_ID %in% setdiff(master.table.raw$QUERY_ORF_ID, lcr.acs),]
#convert to data.frame
master.table.raw<-as.data.frame(master.table.raw)
#assign standard column names
colnames(master.table.raw)[3:ncol(master.table.raw)]<-sapply(colnames(master.table.raw)[3:ncol(master.table.raw)],function(i){return(paste("X",i,sep=""))})

remove<-head(filter(master.table.raw, ID_QUERY_GENOME==3) %>% select(QUERY_ORF_ID,ID_QUERY_GENOME,X3,X4))
remove<-head(master.table.raw)
remove<-head(normalize.scores(raw.df = head(master.table.raw)))

one<-master.table.raw %>% filter(ID_QUERY_GENOME==1) %>% group_by(ID_QUERY_GENOME)
remove<-apply(one[,3:ncol(one)],1,function(i){
  if(i[1]==max(i)){
    return(FALSE)
  }else{
    return(TRUE)
  }
    })
two.remove<-one[remove,]

master.table.normal.signed<-attach.genomeid.header(df = sign.bh.table(bh.table = normalize.scores(raw.df = restrict.minimal.hits(df = master.table.raw))))
master.table.3.6<-select.genomes(df = master.table.normal.signed, g.ids = c(3,6),sep="X")

legend<-load.legend(legend.file = "legend.csv")
gp.3.6<-gen.pca(df = master.table.3.6,legend = legend,g.ids = c(3,6),sep = "X",var.axes = TRUE)
pdf(file = "drosXwolb.pdf",width = 15, height = 15)
gp.3.6
dev.off()

master.table.3.4<-select.genomes(df = master.table.normal.signed, g.ids = c(3,4),sep="X")
pdf(file = "wolbXwolb.pdf",width = 15, height = 15)
gen.pca(df = master.table.3.4,legend = legend,g.ids = c(3,4),sep = "X")
dev.off()

master.table.7.41<-select.genomes(df = master.table.normal.signed, g.ids = c(7,41),sep="X")
pdf(file = "7.41.pdf",width = 15, height = 15)
gen.pca(df = master.table.7.41,legend = legend,g.ids = c(7,41),sep = "X")
dev.off()

master.table.13.8.7.41<-select.genomes(df = master.table.normal.signed, g.ids = c(13,8,7,41),sep="X")
pdf(file = "13.8.7.41.pdf",width = 15, height = 15)
gen.pca(df = master.table.13.8.7.41,legend = legend,g.ids = c(13,8,7,41),sep = "X")
dev.off()

master.table.13.8.7.41.31.23<-select.genomes(df = master.table.normal.signed, g.ids = c(13,8,7,41,31,23),sep="X")
pdf(file = "13.8.7.41.31.23.pdf",width = 15, height = 15)
gen.pca(df = master.table.13.8.7.41.31.23,legend = legend,g.ids = c(13,8,7,41,31,23),sep = "X")
dev.off()

master.table.3.4.1.6<-select.genomes(df = master.table.normal.signed, g.ids = c(3,4,1,6),sep="X")
pdf(file = "wolbXwolbXecoliXdros.pdf",width = 15, height = 15)
gen.pca(df = master.table.3.4.1.6,legend = legend,g.ids = c(3,4,1,6),sep = "X")
dev.off()

master.table.3.4.1.6.14<-select.genomes(df = master.table.normal.signed, g.ids = c(3,4,1,6,14),sep="X")
pdf(file = "wolbXwolbXecoliXdrosXanoph.pdf",width = 15, height = 15)
gen.pca(df = master.table.3.4.1.6.14,legend = legend,g.ids = c(3,4,1,6,14),sep = "X")
dev.off()

master.table.3.4.1.6.14<-select.genomes(df = master.table.normal.signed, g.ids = c(3,4,1,6,14),sep="X")
pdf(file = "wolbXwolbXecoliXdrosXanoph.pdf",width = 15, height = 15)
gen.pca(df = master.table.3.4.1.6.14,legend = legend,g.ids = c(3,4,1,6,14),sep = "X")
dev.off()

old.master.table<-load.hits(hits.file = "master_table5_1.txt")
old.master.table.3.6<-select.genomes(df = old.master.table,g.ids = c(3,6),sep="X")

remove<-getMLE(df = t(select.genomes(df = master.table.3.6,g.ids = 3,sep="X")))

master.table.3.6<-rbind(master.table.3.6,
                        "3XMLE"=getMLE(df = t(select.genomes(df = master.table.3.6,g.ids = 3,sep="X")))$gprime,
                        "6XMLE"=getMLE(df = t(select.genomes(df = master.table.3.6,g.ids = 6,sep="X")))$gprime)

pdf(file = "drosXwolb_MLE.pdf",width = 15, height = 15)
gen.pca(df = master.table.3.6,legend = legend,g.ids = c(3,6),sep = "X", MLE = c(1194:1195))
dev.off()

tmp<-master.table.raw %>% data.frame() %>% restrict.minimal.hits(minhit = 5) %>% normalize.scores() %>% group_by(ID_QUERY_GENOME) %>% split(f=.$ID_QUERY_GENOME)
mles.raw<-lapply(tmp,FUN = function(i){
  message(paste("Calculating MLE for genome",i[1,2]))
  return(getMLE(t(i[,3:ncol(i)])))
})


mle.df<-extract.gprimes(mles.raw)
mle.df<-assign.names(mle.df,legend)
rownames(mle.df)<-colnames(mle.df)
mle.df.norm<-normalize.MLE(mle.df = mle.df)

mle.nj<-nj(X = dist(mle.df))
pdf(file = "tree.nj.pdf", width = 15, height = 15)
plot(mle.nj)
dev.off()

mle.nj<-nj(X = dist(mle.df.norm))
pdf(file = "tree.norm.nj.pdf", width = 15, height = 15)
plot(mle.nj)
dev.off()

pdf(file = "tree.hclust.pdf", width = 15, height = 15)
plot(as.phylo(hclust(dist(mle.df),method = "ward.D2")))
dev.off()

pdf(file = "tree.norm.hclust.pdf", width = 15, height = 15)
plot(as.phylo(hclust(dist(mle.df.norm))))
dev.off()



logreg.mismatch.genes(data = master.table.normal.signed,org1 =3,org2 = 6)


filter(legend,!id_genomes %in% c(2,3))

        

