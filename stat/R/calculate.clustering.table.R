calculate.clustering.table<- function(bh.data, rsquared.cutoff, rsquared.table, legend, processors=4, data.folder=data.folder){
        registerDoMC(processors)
        #for each match table find the genes that are in the core and those, that are outsiders
        core.outsider.list<-foreach(x=1:nrow(rsquared.table), .errorhandling = "stop", .verbose = FALSE)%dopar%{
                match.table<- fread(input = data.folder(rsquared.table$match.file.name[x]), sep = "\t", data.table = FALSE, stringsAsFactors = FALSE)
                #extract i and j from the filename
                ij<- strsplit(x = rsquared.table$match.file.name[x], split = ".", fixed=TRUE)
                ij<- strsplit(x = ij[[1]][1], split = "_", fixed=TRUE)
                i<- as.numeric(ij[[1]][1])
                j<- as.numeric(ij[[1]][2])
                #split by genome ids
                i.columns<- match.table %>% filter(grepl(pattern = paste0("\\b",i,"X"),x = name,perl = TRUE))
                j.columns<- match.table %>% filter(grepl(pattern = paste0("\\b",j,"X"),x = name,perl = TRUE))
                #if the rsquared is too low, then just return a nonsep-filled columns
                if(rsquared.table$rsquared[x]<rsquared.cutoff){
                        i.col=rep("nosep",nrow(i.columns))
                        names(i.col)<- i.columns$name
                        j.col=rep("nosep",nrow(j.columns))
                        names(j.col)<- j.columns$name
                        return(list(i.col=i.col, j.col=j.col, i=i, j=j))
                }else{
                        i.col=sapply(i.columns$match, function(k){
                                if(k){
                                        return("core")
                                }else{
                                        return("outsider")
                                }
                        })
                        names(i.col)<- i.columns$name
                        j.col=sapply(j.columns$match, function(k){
                                if(k){
                                        return("core")
                                }else{
                                        return("outsider")
                                }
                        })
                        names(j.col)<- j.columns$name
                        return(list(i.col=i.col, j.col=j.col, i=i, j=j))
                }
        }
        core.outsider.list<-lapply(legend$id_genomes,function(i){
                sub.list<- lapply(core.outsider.list,function(j){
                        if(j$i==i){
                                return(j)
                        }
                })
                return(sub.list[!sapply(sub.list,is.null)])
        })
        core.outsider.list<- core.outsider.list[-length(core.outsider.list)]
        core.outsider.list<- lapply(core.outsider.list, function(x){
                diagonal.element<- list(i.col = rep("nosep",length(x[[1]]$i.col)),j.col = rep("nosep",length(x[[1]]$i.col)), i=x[[1]]$i, j=x[[1]]$i)
                names(diagonal.element$i.col)<- names(x[[1]]$i.col)
                names(diagonal.element$j.col)<- names(x[[1]]$i.col)
                return(append(x = x,values = list(diagonal.element),0))
        })
        last.genome<- select.genomes(df = bh.data, g.ids = legend$id_genomes[nrow(legend)])
        last.genome.list<-list(list(i.col=rep("nosep", nrow(last.genome)), 
                                    j.col=rep("nosep",nrow(last.genome)), 
                                    i=legend$id_genomes[nrow(legend)], 
                                    j=legend$id_genomes[nrow(legend)]))
        names(last.genome.list[[1]]$i.col)<- rownames(last.genome)
        names(last.genome.list[[1]]$j.col)<- rownames(last.genome)
        core.outsider.list[[length(core.outsider.list)+1]]<- last.genome.list
        return(core.outsider.list)
}