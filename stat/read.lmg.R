if(!file.exists("master.table.core.txt")){
  #generate a matrix big enough to hold the data (essentially the same size as the bh.data), initially filled in with "outsider"
  master.table.core.0.75.cut<-as.data.frame(matrix(nrow=nrow(bh.data),ncol=nrow(legend),factor.core[3]),stringsAsFactors = FALSE)
  #copy the row names from bh.data
  row.names(master.table.core.0.75.cut)<-row.names(bh.data)
  #copy colnames fromt he legend
  colnames(master.table.core.0.75.cut)<-legend$id_genomes
  #now, for each row in the initial bh.table go over the models one by one, load them into ram and see where the row gene is placed 
  #relative to each genome
  registerDoMC(12)
  log.file<-"cluster.log"
  if(file.exists(log.file)){
    file.remove(log.file)
  }
  master.table.core.0.75.cut.list<-foreach(i=1:nrow(legend), .errorhandling = "stop", .verbose = TRUE)%dopar%{
    master.table.core.0.75.cut.i.j<-foreach(j = i:nrow(legend),.verbose = TRUE)%do%{
      #i==j can not be separated, so just return a nonsep list
      sink(log.file, append=TRUE);cat(paste("Doing", i, j, "pair.\n"));sink()
      #load the corresponding model
      if(j!=i){
        load(file = paste(dir,paste(paste(legend$id_genomes[i],legend$id_genomes[j],sep="_"),".rda",sep = ""),sep = "/"))
      }
      sink(log.file, append=TRUE);cat(paste("Loaded", i, j, "pair.\n"));sink()
      #see if the model R squared score if over 0.75, which makes sence in terms of true positive ratio (needs link)
      if(j==i||as.numeric(l.m.g$Rsq)<0.25){
        genome1<-select.genomes(df = bh.data, g.ids = c(legend$id_genomes[i]))
        genome1.is<-rep("nosep", nrow(genome1))
        names(genome1.is)<-rownames(genome1)
        genome2<-select.genomes(df = bh.data, g.ids = c(legend$id_genomes[j]))
        genome2.js<-rep("nosep", nrow(genome2))
        names(genome2.js)<-rownames(genome2)
        sink(log.file, append=TRUE);cat(paste("Returning artificial for", i, j, "pair.\n"));sink()
        return(list(is=genome1.is, js=genome2.js, i=i, j=j))
      }else{
        #for the pair of genomes - restrict the bh.table to only those selected
        data.tab<-select.genomes(df = bh.data, g.ids = c(legend$id_genomes[i],legend$id_genomes[j]))
        #attempt to predict the probabilities for each pair with a precalculated logistic regression
        glm.probs<-predict(l.m.g$fit, data.tab, type = "response")
        #actually split into three groups with a 0.25, 0.25:0.75, 0.75 probability groups
        glm.probs.groups<-sapply(glm.probs,function(t){
          if(t>=0.75){
            return(FALSE)
          }else {
            return(TRUE)
          }
        })
        #append another column, which can be either true of false to indicate the clusters 
        data.tab$cluster<-grepl(pattern = paste("\\b",legend$id_genomes[i],"X",sep = "",collapse = ""),x = rownames(data.tab),perl = TRUE)
        #see if the clustering assignment matches the true clusters, mark as core if yes, otherwise - outsider
        sep.factor.v<-sapply(1:length(glm.probs.groups), function(t){
          if(glm.probs.groups[t]==data.tab$cluster[t]){
            return(factor.core[1])
          }else{ 
            return(factor.core[2])
          }
        })
        names(sep.factor.v)<-names(glm.probs.groups)
        #assign common name identificators like "X1,2..." in accordance with the legend
        sep.factor.v.j<-sep.factor.v[grepl(pattern = paste("\\b",legend$id_genomes[j],"X",sep = "",collapse = ""),x = names(sep.factor.v),perl = TRUE)]
        sep.factor.v.i<-sep.factor.v[grepl(pattern = paste("\\b",legend$id_genomes[i],"X",sep = "",collapse = ""),x = names(sep.factor.v),perl = TRUE)]
        core.col.i<-which(legend$id_genomes==legend$id_genomes[i])
        core.col.j<-which(legend$id_genomes==legend$id_genomes[j])
        #mark i and j cells (because the matrix here is not triangle and needs symmetry)
        master.table.core.0.75.cut.i<-foreach(k=1:length(sep.factor.v.i),.combine = "c")%do%{
          return(master.table.core.0.75.cut[names(sep.factor.v.i)[k],core.col.j]<-sep.factor.v.i[k])
        }
        master.table.core.0.75.cut.j<-foreach(k=1:length(sep.factor.v.j),.combine = "c")%do%{
          return(master.table.core.0.75.cut[names(sep.factor.v.j)[k],core.col.i]<-sep.factor.v.j[k])
        }
        sink(log.file, append=TRUE);cat(paste("Returning real for", i, j, "pair.\n"));sink()
        return(list(is=master.table.core.0.75.cut.i,js=master.table.core.0.75.cut.j,i=core.col.i,j=core.col.j))
      }
    }
    sink(log.file, append=TRUE);cat(paste("Done full row for", i, "\n"));sink()
    return(master.table.core.0.75.cut.i.j)
  }
  #save the list cuz it takes hellalot of time to calculate
  save(file = "master.table.core.0.75.cut.list.rda", master.table.core.0.75.cut.list)
  #so now to concatinate all i into sum matrices
  master.table.core.0.75.cut.list.cbind.is<-lapply(master.table.core.0.75.cut.list, function(x){
    l.is<-do.call(mapply, c(cbind,lapply(x,function(i){
      return(i$is)
    })))
    return(l.is)
  })
  #js are columns in fact
  master.table.core.0.75.cut.list.cbind.js<-lapply(1:nrow(legend), function(x){
    l.js<-lapply(master.table.core.0.75.cut.list[[x]],function(j){
      return(j$js)
    })
    return(l.js)
  })
  master.table.core.0.75.cut.list.cbind.js<-lapply(1:nrow(legend),function(x){
    l.js<-do.call(mapply, c(cbind,lapply(1:(nrow(legend)-length(master.table.core.0.75.cut.list.cbind.js[[x]])+1),function(j){
      return(master.table.core.0.75.cut.list.cbind.js[[j]][[x-j+1]])
    })))
    return(l.js)
  })
  #now concatenate all into one matrix
  master.table.core.0.75.cut.list.cbind<-lapply(1:nrow(legend),function(i){
    dt<-as.data.frame(t(rbind(master.table.core.0.75.cut.list.cbind.is[[i]],
                              master.table.core.0.75.cut.list.cbind.js[[i]])))
    return(dt)
  })
  master.table.core.0.75.cut<- as.data.frame(rbind_all(master.table.core.0.75.cut.list.cbind))
  head(master.table.core.0.75.cut)
  rownames(master.table.core.0.75.cut)<- unlist(sapply(master.table.core.0.75.cut.list.cbind,rownames))
  colnames(master.table.core.0.75.cut)<-sapply(legend$id_genomes,function(i){return(paste0("X",i))})
  #the last column is not needed
  master.table.core.0.75.cut<-master.table.core.0.75.cut[,-ncol(master.table.core.0.75.cut)]
  head(master.table.core.0.75.cut)
  #finally save out the resulting table
  write.table(x=master.table.core.0.75.cut, file="master.table.core.0.75.cut.txt", sep="\t", quote = FALSE)
}