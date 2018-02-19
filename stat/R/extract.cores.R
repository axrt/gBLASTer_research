extract.cores<- function(bh.data, core.outsider.df, processors=6, core.number.cutoff=-1, log.file="cluster.log"){
  cluster <- makeCluster(processors, outfile = log.file)
  registerDoParallel(cluster)
  #first summarize the "core.outsider.df""
  core.outsider.df.summary<- as.data.frame(t(parApply(cl=cluster, X = core.outsider.df, MARGIN = 1, FUN = function(i){
    sum.core<-length(i[i=="core"])
    sum.outsider<-length(i[i=="outsider"])
    sum.nosep<-length(i[i=="nosep"])
    return(c(core=sum.core, outsider=sum.outsider, nonsep=sum.nosep))
  })))
  #do not forget to cluse the cluster
  stopCluster(cluster)
  core.outsider.df.summary<- core.outsider.df.summary %>% mutate(names=rownames(core.outsider.df.summary)) %>%
    filter(outsider==0)
  if(core.number.cutoff>0){
    core.outsider.df.summary<-  core.outsider.df.summary %>% filter(core>=core.number.cutoff)
  }
  return(core.outsider.df.summary$names)
}