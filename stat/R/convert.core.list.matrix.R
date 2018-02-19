convert.core.list.matrix<- function(clustering.list, legend){
        #so now to concatinate all i into sum matrices
        clustering.list.cbind.is<-lapply(clustering.list, function(x){
                l.is<-do.call(mapply, c(cbind,lapply(x,function(i){
                        return(i$i.col)
                })))
                return(l.is)
        })
        #js are columns in fact,so
        clustering.list.cbind.js<-lapply(1:nrow(legend), function(x){
                l.js<-lapply(clustering.list[[x]],function(j){
                        return(j$j.col)
                })
                return(l.js)
        })
        clustering.list.cbind.js<-lapply(1:nrow(legend),function(x){
                l.js<-do.call(mapply, c(cbind,lapply(1:(nrow(legend)-length(clustering.list.cbind.js[[x]])+1),function(j){
                        return(clustering.list.cbind.js[[j]][[x-j+1]])
                })))
                return(l.js)
        })
        #now concatenate all into one matrix
        clustering.list.cbind<-lapply(1:nrow(legend),function(i){
                dt<-as.data.frame(t(rbind(clustering.list.cbind.js[[i]],
                                          clustering.list.cbind.is[[i]])))
                dt<-dt[,-i]
                colnames(dt)<- sapply(legend$id_genomes, function(x){
                  return(paste0("X",x))
                })
                return(dt)
        })
        print(class(clustering.list.cbind)) 
        full.data.frame<- as.data.frame(rbind_all(clustering.list.cbind))
        rownames(full.data.frame)<- unlist(sapply(clustering.list.cbind,rownames))
      
        return(full.data.frame)
}