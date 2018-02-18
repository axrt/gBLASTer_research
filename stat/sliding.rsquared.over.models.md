# Sliding Rsquared Over Models

Lets' see how many logistic models do we have:


```r
if(!require("doParallel")){
  install.packages("doParallel")
  library("doParallel")
}
```

```
## Loading required package: doParallel
## Loading required package: foreach
## Loading required package: iterators
## Loading required package: parallel
```

```r
if(!require("ggplot2")){
  install.packages("ggplot2")
  library("ggplot2")
}
```

```
## Loading required package: ggplot2
```

```r
if(!require("scales")){
  install.packages("scales")
  library("scales")
}
```

```
## Loading required package: scales
```

```r
if(!require("stringr")){
  install.packages("stringr")
  library("stringr")
}
```

```
## Loading required package: stringr
```

```r
if(!require("doMC")){
  install.packages("doMC")
  library("doMC")
}
```

```
## Loading required package: doMC
```

```r
if(!require("foreach")){
  install.packages("foreach")
  library("foreach")
}
if(!require("data.table")){
  install.packages("data.table")
  library("data.table")
}
```

```
## Loading required package: data.table
```

```r
source("R/rphylip.lib.R")
```

```
## Loading required package: Rphylip
## Loading required package: ape
```

```r
source("R/httr.lib.R")
```

```
## Loading required package: httr
```

```r
library("gbra")
```

```
## Loading required package: plyr
## Loading required package: dplyr
## 
## Attaching package: 'dplyr'
## 
## The following objects are masked from 'package:plyr':
## 
##     arrange, count, desc, failwith, id, mutate, rename, summarise,
##     summarize
## 
## The following objects are masked from 'package:data.table':
## 
##     between, last
## 
## The following objects are masked from 'package:stats':
## 
##     filter, lag
## 
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
## 
## Loading required package: tidyr
## Loading required package: Rcpp
## Loading required package: boot
```


```r
#load the helper
source("helper.R")
source("R/helper.R")
#assign self a directory
self.dir<- "sliding.rsquared.over.models.output"
if(!file.exists(self.dir)){
  dir.create(self.dir)
}
self.dir.folder<- function(file=""){
  return(paste(self.dir, file, sep="/"))
}
```


```r
bh.data <- read.table(file = "bh.data.minhit_4.txt") #data
legend<-load.legend(legend.file = "/home/alext/Developer/gbra/inst/extradata/legend.csv",sep = "\t") #legend
```

##Reading models

First let's read all models and see which Rsquareds do they have


```r
dir<-"data" #data folder
#read all models and write the rsuqareds in one table
model.names<- unlist(lapply(1:nrow(legend),function(i){
  return(lapply(i:nrow(legend),function(j){
    if(i!=j){
          return(paste0(paste(legend$id_genomes[i],legend$id_genomes[j],sep="_"),".rda"))
    }
  }))
}))
model.df.file<- "model.df.txt"
if(!file.exists(model.df.file)){
  cluster<- makeCluster(4)
  registerDoParallel(cl=cluster, cores = 4)
  clusterExport(cluster, c("model.names","dir"))
  model.rsquared<- parSapply(cl = cluster, model.names, function(x){
    load(data.folder(x))
    return(l.m.g$Rsq)
  })
  stopCluster(cluster)
  #now get put the table together
  model.df<- data.frame(file.name=model.names, rsquared=model.rsquared, stringsAsFactors = FALSE)
  #save the table
  write.table(x=model.df, file=model.df.file, quote = FALSE, sep="\t")
}
model.df<-read.table(file = model.df.file, sep="\t", stringsAsFactors = FALSE, header = TRUE)
```

Let's do a quick overview of the Rsquared distribution.


```r
#plot the distribution of rsquareds
model.df %>% arrange(rsquared) %>% ggplot(data=., mapping=aes(x=rescale_max(as.numeric(factor(file.name, levels=file.name))), y= rsquared)) + 
  geom_line(color="red") + xlab("Model number (scaled)") + ylab("R squared") + theme_bw()
```

![](sliding.rsquared.over.models_files/figure-html/plot rsquared-1.png) 

##Predicting models

Now, the idea is to read all models, predict the selection of two genomes and save the predicted clusters.


```r
#read models
registerDoMC(4)
model.predictions<- foreach(x=model.names, .combine = "c", .packages = c("gbra"))%dopar%{
  prediction.file<-data.folder(str_replace_all(string = x, pattern = "rda", replacement = "txt"))
  #if the output does not exist
  if(!file.exists(prediction.file)){
    load(data.folder(x))
    #get the model id pair (still separated by "_" for the time)
    ij<- strsplit(x = x[1],split = ".",fixed = TRUE)[[1]][1]
    #now split into individual ids
    ij<-as.numeric(unlist(strsplit(x = ij, split="_", fixed=TRUE)))
    #and assign to variables for convinience
    i<-ij[1]
    j<-ij[2]
    #for the pair of genomes - restrict the bh.table to only those selected (select genomes that correspond to the ids in other words)
    data.tab<-select.genomes(df = bh.data, g.ids = c(i,j))
    #attempt to predict the probabilities for each pair with a precalculated logistic regression
    glm.probs<-predict(l.m.g$fit, data.tab, type = "response")
    names(glm.probs)<- rownames(data.tab)
    #combine the predicitons into one dataframe to write out
    glm.probs<-cbind(name=names(glm.probs), probs=glm.probs) %>% data.frame(stringsAsFactors = FALSE)
    #write out actually
    write.table(x = glm.probs, file=prediction.file, quote=FALSE, sep="\t")
  }
  return(prediction.file)
}
```

Now, go over all prediction tables, see which genes match the true clusters, and which do not, save as `.match` files (tables).


```r
#read models
registerDoMC(4)
model.match<- foreach(x=model.predictions, .combine = "c", .packages = c("gbra"))%dopar%{
  match.file<-paste0(x,".match")
  #see if the file exists already
  if(!file.exists(match.file)){
    #load the predicitons table
    prediction.data<- fread(input = x, sep = "\t", data.table = FALSE, stringsAsFactors = FALSE)[,2:3]
    #make sure the colnames are there (not nessessary actually, probably should remove)
    colnames(prediction.data)<- c("name","probs")
    #get the model id pair (still separated by "_" for the time)
    ij<- strsplit(x = str_replace_all(string = x,pattern = paste0(dir,"/"),replacement = ""),split = ".",fixed = TRUE)[[1]][1]
    #now split into individual ids
    ij<-as.numeric(unlist(strsplit(x = ij, split="_", fixed=TRUE)))
    #and assign to variables for convinience
    i<-ij[1]
    j<-ij[2]
    #assign clusters in accordance with the ids, FALSE is likely to go first
    prediction.data$cluster<- grepl(pattern = paste("\\b",i,"X",sep = "",collapse = ""), x = prediction.data$name, perl = TRUE)
    #for each gene see if the probability favors it to stay in its own cluster, or sends it to the other one
    prediction.data$prediction<- sapply(prediction.data$probs,function(t){
      return(t<0.5)
    })
    #now see if the predictions match the actual clustering
    prediction.data$match<- sapply(1:nrow(prediction.data),function(k){
      return(prediction.data$cluster[k]==prediction.data$prediction[k])
    })
    #finally write out the updated table file, with match column in it this time
    write.table(x = prediction.data, file=match.file, quote=FALSE, sep="\t", row.names = FALSE)
  }
  #return the filename for the future reuse
  return(match.file)
}
```

Let us do a quick summary over the files and try to plot the data to see how the number of mismatches depends on rsquared value.


```r
#first we gonna need a function that returns a number of mismatches for a given file
get.match.percent<- function(match.file){#take care, the function here assumes that the file is given as a full path, not just the filename
  #first - load the file
  matches<- fread(input = match.file, sep = "\t", data.table = FALSE, stringsAsFactors = FALSE)
  #summarize the table to see how many matches and mismatches do we get
  percent.match<- matches %>% group_by(match) %>% summarize(mismatch=n())
  #finally return the ratio
  trues<- as.numeric(percent.match[percent.match$match==TRUE,2])
  percent.match<- trues/nrow(matches)
  return(percent.match)
}
# a good idea would be to go over the models and append the number of matches to the list of models
# so let's mark the input files in the model list firs
model.df$match.file.name<- str_replace_all(string = model.df$file.name, pattern = "rda", replacement = "txt.match")
#now we can go with a regular foreach over the models and get the percent
registerDoMC(4)
model.df$match.percent<- foreach(x=model.df$match.file.name, .combine = "c", .errorhandling = "stop")%dopar%{
  percent.match<- get.match.percent(data.folder(x))
  return(percent.match)
}
#just a quick check, because some of the values are NA, as there may be a case when all matches are FALSE and then it will be NA/nrow(matches),
if(!assertthat::are_equal(any(is.na(model.df$match.percent)),FALSE)){
  model.df$match.percent[is.na(model.df$match.percent)]<-1
}
#now time to plot
gp<- ggplot(data=model.df, mapping = aes(x=match.percent, y=rsquared)) + geom_point(color="red") + geom_smooth() + theme_bw()
gp
```

```
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
```

![](sliding.rsquared.over.models_files/figure-html/exploratory plot for rsquared - mismatches-1.png) 

I think the trend looks good.


```r
#here we need a function that takes rsquared and the master.df table, and returns a core-outsider-nonsep table for all genes
calculate.clustering.table<- function(bh.data, rsquared.cutoff, rsquared.table, legend, processors=4){
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
#let's test it  out for rsquared 1
test.list<- calculate.clustering.table(bh.data = bh.data, rsquared.cutoff = .99, rsquared.table = model.df, legend=legend, processors = 4) #cuz there are in fact no 1 rsquareds, but some are so close, that they are reflected as such
length(test.list[[1]])
```

```
## [1] 45
```

```r
length(test.list[[2]])
```

```
## [1] 44
```

```r
table(test.list[[1]][[1]]$i.col)
```

```
## 
## nosep 
##  1927
```

```r
test.list[[1]][[22]]$j
```

```
## [1] 22
```

```r
table(test.list[[1]][[22]]$i.col)
```

```
## 
## nosep 
##  1927
```

Seems to work and is pretty fast.. Now we need a function to convert this list into a proper matrix.


```r
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
    dt<-as.data.frame(t(rbind(clustering.list.cbind.is[[i]],
                              clustering.list.cbind.js[[i]])))
    return(dt)
  })
  full.data.frame<- as.data.frame(rbind_all(clustering.list.cbind))
  rownames(full.data.frame)<- unlist(sapply(clustering.list.cbind,rownames))
  colnames(full.data.frame)<- sapply(legend$id_genomes, function(x){
    return(paste0("X",x))
  })
  return(full.data.frame)
}
#let's check it out
test.data.frame<- convert.core.list.matrix(clustering.list = test.list, legend = legend)
head(test.data.frame)
```

```
##         X1    X2    X3    X4    X5    X6    X7    X8    X9   X10   X11
## 1X8  nosep nosep nosep nosep nosep nosep nosep nosep nosep nosep nosep
## 1X10 nosep nosep nosep nosep nosep nosep nosep nosep nosep nosep nosep
## 1X17 nosep nosep nosep nosep nosep nosep nosep nosep nosep nosep nosep
## 1X37 nosep nosep nosep nosep nosep nosep nosep nosep nosep nosep nosep
## 1X39 nosep nosep nosep nosep nosep nosep nosep nosep nosep nosep nosep
## 1X54 nosep nosep nosep nosep nosep nosep nosep nosep nosep nosep nosep
##        X12   X13   X14   X15   X16   X17   X18   X19   X20   X21   X22
## 1X8  nosep nosep nosep nosep nosep nosep nosep nosep nosep nosep nosep
## 1X10 nosep nosep nosep nosep nosep nosep nosep nosep nosep nosep nosep
## 1X17 nosep nosep nosep nosep nosep nosep nosep nosep nosep nosep nosep
## 1X37 nosep nosep nosep nosep nosep nosep nosep nosep nosep nosep nosep
## 1X39 nosep nosep nosep nosep nosep nosep nosep nosep nosep nosep nosep
## 1X54 nosep nosep nosep nosep nosep nosep nosep nosep nosep nosep nosep
##       X23   X24   X25   X26   X27   X28   X30   X31   X32   X33   X34
## 1X8  core nosep nosep nosep nosep nosep nosep nosep nosep nosep nosep
## 1X10 core nosep nosep nosep nosep nosep nosep nosep nosep nosep nosep
## 1X17 core nosep nosep nosep nosep nosep nosep nosep nosep nosep nosep
## 1X37 core nosep nosep nosep nosep nosep nosep nosep nosep nosep nosep
## 1X39 core nosep nosep nosep nosep nosep nosep nosep nosep nosep nosep
## 1X54 core nosep nosep nosep nosep nosep nosep nosep nosep nosep nosep
##        X35   X36   X37  X38   X39   X40   X41   X42   X43   X44   X45
## 1X8  nosep nosep nosep core nosep nosep nosep nosep nosep nosep nosep
## 1X10 nosep nosep nosep core nosep nosep nosep nosep nosep nosep nosep
## 1X17 nosep nosep nosep core nosep nosep nosep nosep nosep nosep nosep
## 1X37 nosep nosep nosep core nosep nosep nosep nosep nosep nosep nosep
## 1X39 nosep nosep nosep core nosep nosep nosep nosep nosep nosep nosep
## 1X54 nosep nosep nosep core nosep nosep nosep nosep nosep nosep nosep
##        X46    NA
## 1X8  nosep nosep
## 1X10 nosep nosep
## 1X17 nosep nosep
## 1X37 nosep nosep
## 1X39 nosep nosep
## 1X54 nosep nosep
```

This finally works. And now we just need a function that summarizes the table and extracts only the core genes.


```r
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
#lets check
test.data.frame.names<- extract.cores(bh.data = bh.data, core.outsider.df= test.data.frame)
test.data.frame.names[1:10]
```

```
##  [1] "1X8"  "1X10" "1X17" "1X37" "1X39" "1X54" "1X59" "1X70" "1X82" "1X91"
```

```r
bh.data.cores<-bh.data[test.data.frame.names,]
head(bh.data.cores)
```

```
##      X1        X10       X11       X12       X13       X14       X15
## 1X8   1 0.09707913 0.0000000 0.1170328 0.1585117 0.0000000 0.0000000
## 1X10  1 0.00000000 0.0000000 0.1246878 0.0000000 0.0000000 0.0000000
## 1X17  1 0.24059903 0.0000000 0.2152728 0.2231292 0.0000000 0.0000000
## 1X37  1 0.46309560 0.0000000 0.4823662 0.4377636 0.4515304 0.0000000
## 1X39  1 0.00000000 0.0000000 0.0000000 0.2079324 0.0000000 0.0000000
## 1X54  1 0.56430558 0.5771349 0.5455031 0.5720605 0.5353543 0.4595558
##            X16       X17        X18 X19        X2        X20       X21
## 1X8  0.1387793 0.0000000 0.04752675   0 0.0896809 0.06299841 0.0000000
## 1X10 0.0000000 0.0000000 0.14844449   0 0.1549656 0.00000000 0.0000000
## 1X17 0.0000000 0.1213822 0.00000000   0 0.0000000 0.16243240 0.0000000
## 1X37 0.0000000 0.1299338 0.00000000   0 0.0000000 0.15306429 0.4140899
## 1X39 0.0000000 0.0000000 0.00000000   0 0.0000000 0.00000000 0.2116275
## 1X54 0.5231215 0.4303101 0.45030563   0 0.5290946 0.41211196 0.4359811
##            X22       X23        X24        X25       X26       X27
## 1X8  0.0000000 0.2901256 0.05268397 0.05514616 0.0000000 0.1046985
## 1X10 0.0000000 0.0000000 0.13655844 0.00000000 0.0000000 0.0000000
## 1X17 0.0000000 0.0000000 0.00000000 0.00000000 0.0000000 0.1593715
## 1X37 0.1497620 0.0000000 0.00000000 0.00000000 0.0000000 0.5220083
## 1X39 0.0000000 0.2033016 0.00000000 0.00000000 0.0000000 0.0000000
## 1X54 0.5780336 0.6227969 0.00000000 0.47686301 0.5947443 0.5744544
##             X28         X3       X30       X31       X32       X33
## 1X8  0.09976251 0.06434301 0.0000000 0.1268990 0.1069454 0.1004435
## 1X10 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.1246878
## 1X17 0.27422372 0.00000000 0.0000000 0.1205092 0.2475824 0.0000000
## 1X37 0.11506626 0.00000000 0.0000000 0.3011966 0.5423940 0.0000000
## 1X39 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000
## 1X54 0.66427537 0.62220036 0.4320997 0.5329682 0.4550779 0.5535602
##            X34       X35        X36       X37       X38       X39
## 1X8  0.0000000 0.0000000 0.06994843 0.0000000 0.2901256 0.0000000
## 1X10 0.0000000 0.0000000 0.00000000 0.0000000 0.0000000 0.0000000
## 1X17 0.0000000 0.0000000 0.00000000 0.0000000 0.0000000 0.0000000
## 1X37 0.0000000 0.0000000 0.00000000 0.0000000 0.0000000 0.0000000
## 1X39 0.0000000 0.0000000 0.00000000 0.0000000 0.2033016 0.0000000
## 1X54 0.4488182 0.2733442 0.50581427 0.3968887 0.6227969 0.5669939
##              X4       X40        X41        X42        X43       X44
## 1X8  0.05671777 0.0000000 0.08945389 0.05044296 0.07353403 0.0000000
## 1X10 0.00000000 0.0000000 0.16685166 0.00000000 0.12763235 0.0000000
## 1X17 0.00000000 0.0000000 0.00000000 0.16723917 0.23097416 0.0000000
## 1X37 0.00000000 0.0000000 0.00000000 0.22850281 0.34689997 0.5165044
## 1X39 0.00000000 0.0000000 0.00000000 0.00000000 0.00000000 0.0000000
## 1X54 0.63293797 0.1065239 0.41927037 0.44195415 0.46880593 0.4589592
##            X45       X46        X5        X6        X7        X8        X9
## 1X8  0.1006647 0.0000000 0.1004319 0.0000000 0.0000000 0.1455081 0.0000000
## 1X10 0.0000000 0.1656646 0.0000000 0.0000000 0.1644775 0.1258749 0.0000000
## 1X17 0.2956274 0.0000000 0.2497704 0.0000000 0.0000000 0.1117800 0.0000000
## 1X37 0.0000000 0.0000000 0.0000000 0.4476705 0.0000000 0.1398550 0.2196966
## 1X39 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
## 1X54 0.5076039 0.4264365 0.6257796 0.0000000 0.4186738 0.6117494 0.4186738
```

Well, that takes some time, but it works, which is important). Now we need a function that extracts the mles. The `getMLEs()` function takes a RAW `data.frame`, so we need to cook it first. which is easy in fact.


```r
#extract mles from the dataframe
bh.data.cores.mles<- getMLEs(df = expand.df(bh.data.cores))
```

```
## Calculating MLE for genome 8
## Calculating MLE for genome 38435
## Calculating MLE for genome 68603
## Calculating MLE for genome 74671
## Calculating MLE for genome 81279
## Calculating MLE for genome 139217
## Calculating MLE for genome 247626
## Calculating MLE for genome 953029
## Calculating MLE for genome 1087433
## Calculating MLE for genome 1798916
## Calculating MLE for genome 1864756
## Calculating MLE for genome 1898953
## Calculating MLE for genome 1961502
## Calculating MLE for genome 2233078
## Calculating MLE for genome 4328559
## Calculating MLE for genome 4350340
## Calculating MLE for genome 4381254
## Calculating MLE for genome 6713786
## Calculating MLE for genome 6742191
## Calculating MLE for genome 6753817
## Calculating MLE for genome 6935691
## Calculating MLE for genome 7636309
## Calculating MLE for genome 8069566
## Calculating MLE for genome 8268079
## Calculating MLE for genome 8285214
## Calculating MLE for genome 8298835
## Calculating MLE for genome 8766415
## Calculating MLE for genome 9094481
## Calculating MLE for genome 9105946
## Calculating MLE for genome 10115503
## Calculating MLE for genome 10613611
## Calculating MLE for genome 11252396
## Calculating MLE for genome 11281406
## Calculating MLE for genome 11283907
## Calculating MLE for genome 11287593
## Calculating MLE for genome 11339100
## Calculating MLE for genome 11358995
## Calculating MLE for genome 11557810
## Calculating MLE for genome 11664730
## Calculating MLE for genome 11772215
## Calculating MLE for genome 12124931
## Calculating MLE for genome 12294164
## Calculating MLE for genome 12430564
## Calculating MLE for genome 16244488
## Calculating MLE for genome 16264321
```

```r
bh.data.cores.mles.gprimes<- extract.gprimes(bh.data.cores.mles)
head(bh.data.cores.mles.gprimes)
```

```
##             X1        X10        X11        X12        X13        X14
## X1 1.000000000 0.04713139 0.01513603 0.04802937 0.05904990 0.03639293
## X2 0.212921109 0.05550690 0.01870837 0.05597110 0.06773842 0.04312719
## X3 0.204897255 0.05509918 0.02047130 0.05686221 0.08508530 0.04760928
## X4 0.192152762 0.05061854 0.02009982 0.05146275 0.08003888 0.04539256
## X5 0.139744067 0.03836201 0.01003550 0.04350713 0.05108954 0.04847089
## X6 0.003833656 0.02554572 0.01171938 0.02414203 0.02244535 0.35645594
##            X15         X16         X17         X18         X19          X2
## X1 0.010049000 0.065407093 0.005456173 0.059877879 0.035647920 0.179677359
## X2 0.013351152 0.086918289 0.007244052 0.074957681 0.046561763 1.000000000
## X3 0.016737526 0.062006707 0.011079283 0.048444513 0.026678525 0.169160293
## X4 0.015559192 0.059931988 0.010747056 0.044532439 0.024583145 0.154055065
## X5 0.007198221 0.054697062 0.005453318 0.048000942 0.034080712 0.110409453
## X6 0.005572828 0.003738292 0.059435987 0.001965457 0.001306599 0.005156893
##           X20        X21        X22        X23         X24         X25
## X1 0.01813680 0.02554019 0.01517226 0.06315137 0.046107984 0.047180441
## X2 0.02012269 0.03157624 0.01755038 0.07476817 0.056939496 0.066200720
## X3 0.02241434 0.03535106 0.02281367 0.09114785 0.033671624 0.043014146
## X4 0.02045816 0.03626929 0.02231612 0.08572662 0.030776536 0.039923930
## X5 0.01584942 0.03143186 0.01843766 0.05864740 0.039009719 0.037135986
## X6 0.01756645 0.27755784 0.03713486 0.01682451 0.001748482 0.001961275
##           X26        X27         X28          X3         X30        X31
## X1 0.02385326 0.03717991 0.111534112 0.064081349 0.007459109 0.05675107
## X2 0.02805731 0.04173111 0.106948419 0.068837438 0.008724714 0.06477240
## X3 0.03041508 0.04110969 0.269915427 1.000000000 0.009521486 0.04623021
## X4 0.02728543 0.03896918 0.256290505 0.808282414 0.009967763 0.04391263
## X5 0.02665485 0.04572809 0.092417230 0.060975961 0.003790641 0.06899589
## X6 0.03118041 0.04312706 0.001549274 0.001370133 0.038285458 0.03487767
##           X32         X33          X34          X35         X36
## X1 0.04093845 0.148920463 0.0175343488 0.0197396737 0.118996866
## X2 0.04608843 0.197041411 0.0266159758 0.0312943973 0.148673240
## X3 0.04849457 0.145001067 0.0330085370 0.0378044274 0.131414016
## X4 0.04397139 0.133460915 0.0302142472 0.0343403984 0.120163138
## X5 0.04403582 0.124468984 0.0078501977 0.0090914791 0.118413765
## X6 0.03419909 0.004373952 0.0005407261 0.0005070277 0.003796229
##            X37        X38        X39          X4         X40        X41
## X1 0.008874474 0.06315137 0.02132772 0.062773272 0.001500834 0.02992039
## X2 0.011934066 0.07476817 0.02789823 0.066928184 0.001686256 0.03782501
## X3 0.015369419 0.09114785 0.04244075 0.830035111 0.001827472 0.05261133
## X4 0.013845150 0.08572662 0.04128897 0.999995651 0.002026043 0.04809886
## X5 0.005183036 0.05864740 0.01594879 0.059602334 0.001373541 0.02932006
## X6 0.010078203 0.01682451 0.01500524 0.001415788 0.001771402 0.02408390
##           X42        X43        X44         X45        X46          X5
## X1 0.04378764 0.02468855 0.04011100 0.048668635 0.02648413 0.175069466
## X2 0.05477825 0.02934860 0.04506462 0.062129715 0.03398480 0.161588814
## X3 0.06078434 0.02408986 0.08796710 0.037849149 0.03806326 0.297019648
## X4 0.05930378 0.02299410 0.11015421 0.033870556 0.03641253 0.282555731
## X5 0.04339665 0.01514917 0.05065295 0.046761019 0.02702488 0.999994451
## X6 0.02381645 0.02065398 0.36344932 0.002830105 0.02715466 0.004152214
##             X6         X7         X8         X9
## X1 0.007300909 0.02922502 0.07038928 0.01787457
## X2 0.008465303 0.03744179 0.07864243 0.02260965
## X3 0.008555660 0.04424241 0.10392863 0.02357271
## X4 0.010123276 0.04179428 0.09646876 0.02112783
## X5 0.006383382 0.02765162 0.06217918 0.02135671
## X6 1.000000000 0.02641569 0.02259168 0.05784408
```

Now let's see if these gprimes make any sense. We first convert them into a `dist` object with Euclidean distance, second we will feed the output to `fitch`, third we will visualize the data.


```r
#calculate the distance
bh.data.cores.mles.gprimes.dist<- dist(bh.data.cores.mles.gprimes)
#run fitch with a bunch of default parameters (need to comment on this a little more explicitly)
rf<- Rfitch(D=bh.data.cores.mles.gprimes.dist,path = "/usr/local/bin", model="FM", power=2, negative=TRUE, global=FALSE,random.order=FALSE, root=FALSE, quiet=TRUE)
```

```
## arguments 'show.output.on.console', 'minimized' and 'invisible' are for Windows only
```

```
## 
##   45 Populations
## 
## Fitch-Margoliash method version 3.696
## 
##                   __ __             2
##                   \  \   (Obs - Exp)
## Sum of squares =  /_ /_  ------------
##                                 2
##                    i  j      Obs
## 
## Negative branch lengths not allowed
## 
## 
##   +----------------2         
##   ! 
##   !       +---------------28        
##   !    +-26  
##   !    !  !                 +---4         
##   !    !  +-----------------2 
##   !  +-3                    +---3         
##   !  ! ! 
##   !  ! !       +------------34        
##   !  ! +------32  
##   !  !         +------------33        
##   !  !  
##   !  !    +------------------32        
##   !  !    !  
##   !  !    !  +------------------35        
##   !  !    !  !  
##   !  !    !  !     +-----------------44        
##   !  ! +-30  !     !  
##   !  ! !  !  !  +-42     +---------------24        
##   !  ! !  !  !  !  !  +-22  
##   !  ! !  !  !  !  !  !  +---------------19        
##   !  ! !  !  !  !  +-17  
##   !  ! !  !  !  !     !  +-----------------18        
##   !  ! !  +-33  !     +-16  
##   !  ! !     !  !        !  +---------------25        
##   !  ! !     !  !        +-23  
##   !  ! !     !  !           +---------------16        
##   !  ! !     !  !  
##   !  ! !     !  !                             +37        
##   1-31 !     !  ! +--------------------------35  
##   !  ! !     !  ! !                           +23        
##   !  ! !     !  ! ! 
##   !  ! !     !  ! !                  +-----------40        
##   !  ! !     +-14 !      +----------38  
##   !  ! !        ! !      !           !  +----------45        
##   !  ! !        ! !      !           +-43  
##   !  ! !        ! !   +-11              +-----------7         
##   !  ! !        ! !   !  !  
##   !  ! !        ! !   !  !  +------------------27        
##   !  ! !        ! !   !  !  !  
##   !  ! !        ! !   !  +-25  +-------------------26        
##   !  ! !        ! !   !     !  !  
##   !  ! !        ! !   !     +-24  +---------------------39        
##   !  ! !        ! !   !        !  !  
##   !  ! !        ! !   !        +-37  +------------------22        
##   !  ! !        ! !   !           !  !  
##   !  ! !        ! !   !           +-20 +-----------------9         
##   !  ! !        +-8   !              ! ! 
##   !  ! !          !   !              ! !     +------------------29        
##   !  ! !          !   !              +-7  +-27  
##   !  ! !          ! +-9                !  !  +-----------------17        
##   !  +-4          ! ! !                !  !  
##   !    !          ! ! !                +-15     +---------------21        
##   !    !          ! ! !                   !  +-19  
##   !    !          ! ! !                   !  !  !  +------------43        
##   !    !          ! ! !                   +-12  +-41  
##   !    !          ! ! !                      !     +-------------14        
##   !    !          ! ! !                      !  
##   !    !          ! ! !                      +-----------------6         
##   !    !          ! ! ! 
##   !    !          ! ! !     +------------------38        
##   !    !          ! ! !  +-36  
##   !    !          ! ! !  !  !  +------------------42        
##   !    !          ! ! !  !  +-40  
##   !    !          ! ! !  !     +------------------11        
##   !    !          ! ! +-21  
##   !    !          +-5    !  +-------------------41        
##   !    !            !    !  !  
##   !    !            !    +-39      +----------------36        
##   !    !            !       !  +--34  
##   !    !            !       !  !   +----------------15        
##   !    !            !       +-13  
##   !    !            !          !  +-----------------31        
##   !    !            !          !  !  
##   !    !            !          +-29     +-----------------30        
##   !    !            !             !  +-28  
##   !    !            !             !  !  +------------------20        
##   !    !            !             +-18  
##   !    !            !                !  +----------------12        
##   !    !            !                +-10  
##   !    !            !                   +----------------10        
##   !    !            ! 
##   !    !            ! +------------------13        
##   !    !            +-6 
##   !    !              +------------------8         
##   !    ! 
##   !    +-----------------5         
##   ! 
##   +----------------1         
## 
## 
## remember: this is an unrooted tree!
## 
## Sum of squares =     0.48994
## 
## Average percent standard deviation =     1.57384
## 
## Between        And            Length
## -------        ---            ------
##    1          2                 0.56870
##    1            31              0.02138
##   31             3              0.02377
##    3            26              0.09796
##   26          28                0.53026
##   26             2              0.58877
##    2          4                 0.12284
##    2          3                 0.13705
##    3            32              0.25785
##   32          34                0.43193
##   32          33                0.43556
##   31             4              0.00949
##    4            30              0.01941
##   30          32                0.61033
##   30            33              0.00688
##   33          35                0.61670
##   33            14              0.04625
##   14            42              0.06010
##   42          44                0.60405
##   42            17              0.01772
##   17            22              0.08252
##   22          24                0.53397
##   22          19                0.54124
##   17            16              0.01332
##   16          18                0.58591
##   16            23              0.07509
##   23          25                0.53231
##   23          16                0.51417
##   14             8              0.06954
##    8            35              0.91382
##   35          37                0.00025
##   35          23                0.00021
##    8             5              0.01132
##    5             9              0.01297
##    9            11              0.00848
##   11            38              0.37302
##   38          40                0.40986
##   38            43              0.03416
##   43          45                0.36485
##   43          7                 0.38096
##   11            25              0.01737
##   25          27                0.63411
##   25            24              0.00763
##   24          26                0.64848
##   24            37              0.00739
##   37          39                0.72639
##   37            20              0.04916
##   20          22                0.63436
##   20             7              0.01624
##    7          9                 0.60362
##    7            15              0.02294
##   15            27              0.01959
##   27          29                0.62618
##   27          17                0.59834
##   15            12              0.07014
##   12            19              0.02366
##   19          21                0.52194
##   19            41              0.06311
##   41          43                0.42544
##   41          14                0.44553
##   12          6                 0.60764
##    9            21              0.00302
##   21            36              0.03471
##   36          38                0.63496
##   36            40              0.01205
##   40          42                0.61484
##   40          11                0.61954
##   21            39              0.00614
##   39          41                0.63735
##   39            13              0.00944
##   13            34              0.11463
##   34          36                0.56249
##   34          15                0.56073
##   13            29              0.03011
##   29          31                0.60575
##   29            18              0.03496
##   18            28              0.00793
##   28          30                0.57363
##   28          20                0.61728
##   18            10              0.03080
##   10          12                0.55189
##   10          10                0.54945
##    5             6              0.01870
##    6          13                0.62379
##    6          8                 0.61399
##    4          5                 0.59197
##    1          1                 0.57093
```

```r
rf.names<- rename.tree(tree = rf, legend = legend)
write.tree(phy = rf.names, file = self.dir.folder("test.newick"))
```

Unfortunately, there is no easy way to plot the `fitch`-generated object, so i guess we are going to use iTOL for that.


```r
plot.itol<- function(newick.file, project.name="Default", output.file){
  #request
  req<-POST(url = "http://itol.embl.de/batch_uploader.cgi",
            encode = "multipart", body=list(
              treeFile=upload_file(newick.file),#the tree itself in newick format
              treeFormat="newick", #this says that the tree is in newick and is best to indicate explicitly
              projectName="Default" #just a project name, mandatory though
            ))
  id<- str_replace_all(string = str_trim(toString(req)),pattern = "SUCCESS: ", replacement = "")#this is simply to extract the ID from the server responce
  print(id)
  #download
  pl<-POST(url = "http://itol.embl.de/batch_downloader.cgi",
           encode = "multipart", body=list(
             tree= id,
             format="pdf",#we want it in pdf
             displayMode="unrooted",#we want it in circular and unrooted
             resolution=300,#print resolution, makes no sence for a vector pdf file, but there might ne smth like fonts for internal nodes, that get rasterized
             fontSize=48,#this one is the highest number of the font size, otherwise the letters start to overlap
             lineWidth=1,#this is the line stoke, enough to reflect nicely
             hideRanges=0,#we do not have any ranges here, so it is good idea to switch this off
             scaleFactor=0.8#this is the minimum scaling factor, if made smaller, the words begin to overlap (scaling here is the width of the plot if looked at horyzontally)
           ))
  #extract the binary object from the responce and save it as pdf (which it is).
  writeBin(pl$content, con = output.file)
  return(output.file)
}
#let's try if tthe above works as expected
plot.itol(newick.file = self.dir.folder("test.newick"), output.file = self.dir.folder("test.newick.pdf"))
```

```
## [1] "1291711501382474114412892860"
```

```
## [1] "sliding.rsquared.over.models.output/test.newick.pdf"
```
<embed src="sliding.rsquared.over.models.output/test.newick.pdf" width="900" height="500" type='application/pdf'/>

Eeeeehm.. unfortunately, iTOL is not capable of producinig a clear image either, so I guess the only option is still [Archaeopteryx](https://aptxevo.wordpress.com/). So, let's keep it simple for now and save only the newicks. Let's check out what the trees will look like when we go over (0,1) rsquared (yes those are curve brackets as I do not think there may exist a true 0 or a true 1 among the models) with a step of 0.1.


```r
#create a sequence of rsquared values
rsquareds<- seq(from=0.1, to=0.99, by=0.98/10)
processors<- 6
slide.over.rsquared<- function(rsquareds, cores, legend, bh.data, model.df){
  output.names<- sapply(rsquareds,function(r){
    output.name<- self.dir.folder(paste0(r,".rsq.cut.newick"))
    print(paste("Doing: ",output.name))
    if(!file.exists(output.name)){
      print("Clustering list")
      clustering.list<- calculate.clustering.table(bh.data = bh.data, rsquared.cutoff = r, rsquared.table = model.df, legend = legend, processors = processors)
      print("Merging clustering list")
      clustering.table<- convert.core.list.matrix(clustering.list = clustering.list, legend = legend)
      print("Extracting cores")
      bh.data.cores<- bh.data[extract.cores(bh.data = bh.data, core.outsider.df = clustering.table, processors = processors),]
      print("Extracting MLEs")
      bh.data.cores.mles<- getMLEs(df = expand.df(bh.data.cores))
      print("Extracting gprimes")
      bh.data.cores.mles.gprimes<- extract.gprimes(mle.list = bh.data.cores.mles)
      print("Calculating distance")
      bh.data.cores.mles.gprimes.dist<- dist(bh.data.cores.mles.gprimes)
      print("Running fitch")
      rf<- Rfitch(D=bh.data.cores.mles.gprimes.dist,path = "/usr/local/bin", model="FM", power=2, negative=TRUE, global=FALSE,random.order=FALSE, root=FALSE, quiet=TRUE)
      print("Renaming edges")
      rf.names<- rename.tree(tree = rf, legend = legend)
      print("Saving newick data")
      write.tree(phy = rf.names, file = output.name)
    }
    print(paste("Done: ",output.name))
    return(output.name)
  })
  return(output.names)
}
output.names<- slide.over.rsquared(rsquareds = rsquareds, cores = cores, legend=legend, bh.data=bh.data, model.df = model.df)
```

```
## [1] "Doing:  sliding.rsquared.over.models.output/0.1.rsq.cut.newick"
## [1] "Done:  sliding.rsquared.over.models.output/0.1.rsq.cut.newick"
## [1] "Doing:  sliding.rsquared.over.models.output/0.198.rsq.cut.newick"
## [1] "Done:  sliding.rsquared.over.models.output/0.198.rsq.cut.newick"
## [1] "Doing:  sliding.rsquared.over.models.output/0.296.rsq.cut.newick"
## [1] "Done:  sliding.rsquared.over.models.output/0.296.rsq.cut.newick"
## [1] "Doing:  sliding.rsquared.over.models.output/0.394.rsq.cut.newick"
## [1] "Done:  sliding.rsquared.over.models.output/0.394.rsq.cut.newick"
## [1] "Doing:  sliding.rsquared.over.models.output/0.492.rsq.cut.newick"
## [1] "Done:  sliding.rsquared.over.models.output/0.492.rsq.cut.newick"
## [1] "Doing:  sliding.rsquared.over.models.output/0.59.rsq.cut.newick"
## [1] "Done:  sliding.rsquared.over.models.output/0.59.rsq.cut.newick"
## [1] "Doing:  sliding.rsquared.over.models.output/0.688.rsq.cut.newick"
## [1] "Done:  sliding.rsquared.over.models.output/0.688.rsq.cut.newick"
## [1] "Doing:  sliding.rsquared.over.models.output/0.786.rsq.cut.newick"
## [1] "Done:  sliding.rsquared.over.models.output/0.786.rsq.cut.newick"
## [1] "Doing:  sliding.rsquared.over.models.output/0.884.rsq.cut.newick"
## [1] "Done:  sliding.rsquared.over.models.output/0.884.rsq.cut.newick"
## [1] "Doing:  sliding.rsquared.over.models.output/0.982.rsq.cut.newick"
## [1] "Done:  sliding.rsquared.over.models.output/0.982.rsq.cut.newick"
```

Ok, we are done, now let's see what we got in the [folder](sliding.rsquared.over.models.output/).
So, apparently, the best tree in our perception we get when we use **Rsquared in range [0.7,0.78]** as there we see a clean separation between domains. The thing we would like to do now is throw out the *Ichthyophthirius multifiliis* as we had screwd up with his genetic code, and never did a proper conversion for it.


```r
output.name<- self.dir.folder(paste0("0.78.no.ich.mult.no.ciona",".rsq.cut.newick"))
print(paste("Doing: ",output.name))
```

```
## [1] "Doing:  sliding.rsquared.over.models.output/0.78.no.ich.mult.no.ciona.rsq.cut.newick"
```

```r
if(!file.exists(output.name)){
  print("Clustering list")
  clustering.list<- calculate.clustering.table(bh.data = bh.data, rsquared.cutoff = 0.78, rsquared.table = model.df, legend = legend, processors = processors)
  print("Merging clustering list")
  clustering.table<- convert.core.list.matrix(clustering.list = clustering.list, legend = legend)
  print("Extracting cores")
  bh.data.cores<- bh.data[extract.cores(bh.data = bh.data, core.outsider.df = clustering.table, processors = processors),]
  #remove the Ichthyophthirius multifiliis thing
  bh.data.cores<- bh.data.cores[!grepl(x = rownames(bh.data.cores), pattern = "\\b40X",fixed = FALSE),-which(colnames(bh.data.cores)=="X40")]
  bh.data.cores<- bh.data.cores[!grepl(x = rownames(bh.data.cores), pattern = "\\b22X",fixed = FALSE),-which(colnames(bh.data.cores)=="X22")]
  print("Extracting MLEs")
  bh.data.cores.mles<- getMLEs(df = expand.df(bh.data.cores))
  print("Extracting gprimes")
  bh.data.cores.mles.gprimes<- extract.gprimes(mle.list = bh.data.cores.mles)
  print("Calculating distance")
  bh.data.cores.mles.gprimes.dist<- dist(bh.data.cores.mles.gprimes)
  print("Running fitch")
  rf<- Rfitch(D=bh.data.cores.mles.gprimes.dist,path = "/usr/local/bin", model="FM", power=2, negative=TRUE, global=FALSE,random.order=FALSE, root=FALSE, quiet=TRUE)
  print("Renaming edges")
  rf.names<- rename.tree(tree = rf, legend = legend)
  print("Saving newick data")
  write.tree(phy = rf.names, file = output.name)
  plot.itol(newick.file = output.name, output.file = paste0(output.name,".pdf"))
}
print(paste("Done: ",output.name))
```

```
## [1] "Done:  sliding.rsquared.over.models.output/0.78.no.ich.mult.no.ciona.rsq.cut.newick"
```

##Tightening cores

At this point is seems like I am going to do this a lot, so it makes sense to create a function out of the above boilerplate code. But before all these components go into the `gbra` package, let's see if all works well:

Below we are goingt ot tighten the cores, which means that we are planning ot put a rigid cutoff on the number of times a gene has to be in the core to be considered for future tree prediction. Imgaine a situation, when a gene A has `core, nosep, nosep, nosep`. Technically, that allows it to be a core gene. However, it is drammatically different from the following case: `core, core, core, core`, where the gene is clearly a very specific point in the genome array.

let's use Fibonacci numbers, starting with 3 up to 34.


```r
#here goes a function to go over the whole process in one call
build.tree<- function(bh.data, core.number.cutoff=-1, rsquared.cutoff, rsquared.table, legend, processors, output.name){
  print(paste("Doing: ",output.name))
  if(!file.exists(output.name)){
    print("Clustering list")
    clustering.list<- calculate.clustering.table(bh.data = bh.data, rsquared.cutoff = rsquared.cutoff, rsquared.table = model.df, legend = legend, processors = processors)
    print("Merging clustering list")
    clustering.table<- convert.core.list.matrix(clustering.list = clustering.list, legend = legend)
    print("Extracting cores")
    bh.data.cores<- bh.data[extract.cores(bh.data = bh.data, core.number.cutoff=core.number.cutoff, core.outsider.df = clustering.table, processors = processors),]
    print("Extracting MLEs")
    bh.data.cores.mles<- getMLEs(df = expand.df(bh.data.cores))
    print("Extracting gprimes")
    bh.data.cores.mles.gprimes<- extract.gprimes(mle.list = bh.data.cores.mles)
    print("Calculating distance")
    bh.data.cores.mles.gprimes.dist<- dist(bh.data.cores.mles.gprimes)
    print("Running fitch")
    rf<- Rfitch(D=bh.data.cores.mles.gprimes.dist,path = "/usr/local/bin", model="FM", power=2, negative=TRUE, global=FALSE,random.order=FALSE, root=FALSE, quiet=TRUE)
    print("Renaming edges")
    rf.names<- rename.tree(tree = rf, legend = legend)
    print("Saving newick data")
    write.tree(phy = rf.names, file = output.name)
    plot.itol(newick.file = output.name, output.file = paste0(output.name,".pdf"))
  }
  print(paste("Done: ",output.name))
  return(output.name)
}

#let's generate the finonacci sequence
#was lasy, so copypasted from R-bloggers http://www.r-bloggers.com/example-7-1-create-a-fibonacci-sequence/
len <- 9
fibvals <- numeric(len)
fibvals[1] <- 1
fibvals[2] <- 1
for (i in 3:len) { 
   fibvals[i] <- fibvals[i-1]+fibvals[i-2]
} 
fibvals<- fibvals[2:length(fibvals)]
fibvals
```

```
## [1]  1  2  3  5  8 13 21 34
```

```r
#now let's run the function these 5 times..
rsquared.cutoff<-0.73
tree.newick.files<- sapply(fibvals, function(f){
  output.name<- self.dir.folder(paste0("r.",rsquared.cutoff,".corecut.",f,".newick"))
  if(!file.exists(output.name)){
    return(build.tree(
    bh.data = bh.data, core.number.cutoff = f, 
    rsquared.cutoff = rsquared.cutoff, rsquared.table = model.df, 
    processors = 6, legend = legend, output.name = output.name))
  }else{
      return(output.name)
  }
})
```

```
## [1] "Doing:  sliding.rsquared.over.models.output/r.0.73.corecut.1.newick"
## [1] "Clustering list"
## [1] "Merging clustering list"
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## [1] "Extracting cores"
## [1] "Extracting MLEs"
```

```
## Calculating MLE for genome 8
## Calculating MLE for genome 38435
## Calculating MLE for genome 68603
## Calculating MLE for genome 74671
## Calculating MLE for genome 81279
## Calculating MLE for genome 139217
## Calculating MLE for genome 247626
## Calculating MLE for genome 953029
## Calculating MLE for genome 1087433
## Calculating MLE for genome 1798918
## Calculating MLE for genome 1864829
## Calculating MLE for genome 1898956
## Calculating MLE for genome 1961502
## Calculating MLE for genome 2233078
## Calculating MLE for genome 4328604
## Calculating MLE for genome 4350340
## Calculating MLE for genome 4381525
## Calculating MLE for genome 6713786
## Calculating MLE for genome 6742235
## Calculating MLE for genome 6753817
## Calculating MLE for genome 6935691
## Calculating MLE for genome 7636370
## Calculating MLE for genome 8069566
## Calculating MLE for genome 8268079
## Calculating MLE for genome 8285225
## Calculating MLE for genome 8298835
## Calculating MLE for genome 8766415
## Calculating MLE for genome 9094481
## Calculating MLE for genome 9105946
## Calculating MLE for genome 10115503
## Calculating MLE for genome 10613611
## Calculating MLE for genome 11252396
## Calculating MLE for genome 11281406
## Calculating MLE for genome 11283907
## Calculating MLE for genome 11287593
## Calculating MLE for genome 11339100
## Calculating MLE for genome 11358995
## Calculating MLE for genome 11557810
## Calculating MLE for genome 11664862
## Calculating MLE for genome 11772215
## Calculating MLE for genome 12124931
## Calculating MLE for genome 12294189
## Calculating MLE for genome 12430564
## Calculating MLE for genome 16244488
## Calculating MLE for genome 16264321
```

```
## [1] "Extracting gprimes"
## [1] "Calculating distance"
## [1] "Running fitch"
```

```
## arguments 'show.output.on.console', 'minimized' and 'invisible' are for Windows only
```

```
## 
##   45 Populations
## 
## Fitch-Margoliash method version 3.696
## 
##                   __ __             2
##                   \  \   (Obs - Exp)
## Sum of squares =  /_ /_  ------------
##                                 2
##                    i  j      Obs
## 
## Negative branch lengths not allowed
## 
## 
##   +----------------2         
##   ! 
##   !              +-----------34        
##   !    +--------32  
##   !    !         +-----------33        
##   ! +-31  
##   ! !  !   +---------------28        
##   ! !  +--26  
##   ! !      !                 +--4         
##   ! !      +-----------------2 
##   ! !                        +---3         
##   ! ! 
##   ! ! +----------------5         
##   1-3 ! 
##   ! ! !     +-----------------35        
##   ! ! !  +-33  
##   ! ! !  !  +-----------------32        
##   ! ! !  !  
##   ! ! !  !     +-----------------44        
##   ! ! !  !     !  
##   ! ! !  !  +-42     +--------------24        
##   ! +-4  !  !  !  +-22  
##   !   !  !  !  !  !  +---------------19        
##   !   !  !  !  +-17  
##   !   !  !  !     !  +----------------18        
##   !   !  !  !     +-16  
##   !   !  !  !        !  +--------------25        
##   !   !  !  !        +-23  
##   !   !  !  !           +-------------16        
##   !   !  !  !  
##   !   +-30  !           +------------------27        
##   !      !  !           !  
##   !      !  !        +-25  +-------------------26        
##   !      !  !        !  !  !  
##   !      !  !        !  +-24  +-----------------------39        
##   !      !  !        !     !  !  
##   !      !  !        !     +-37  +------------------22        
##   !      !  !        !        !  !  
##   !      !  !        !        +-20 +-----------------9         
##   !      !  !        !           ! ! 
##   !      !  !        !           ! !     +------------------29        
##   !      !  !        !           +-7  +-27  
##   !      !  !        !             !  !  +-----------------17        
##   !      !  !        !             !  !  
##   !      !  !     +-13             +-15     +---------------21        
##   !      !  !     !  !                !  +-19  
##   !      +-14     !  !                !  !  !  +------------43        
##   !         !     !  !                +-12  +-41  
##   !         !     !  !                   !     +-------------14        
##   !         !     !  !                   !  
##   !         !     !  !                   +------------------6         
##   !         !     !  !  
##   !         !     !  !     +------------------38        
##   !         !     !  !  +-36  
##   !         !     !  !  !  !  +-----------------42        
##   !         !     !  !  !  +-40  
##   !         !     !  !  !     +------------------11        
##   !         !     !  +-11  
##   !         !     !     !  +------------------41        
##   !         !     !     !  !  
##   !         !  +-21     +-39      +---------------36        
##   !         !  !  !        ! +---34  
##   !         !  !  !        ! !    +---------------15        
##   !         !  !  !        +-9 
##   !         !  !  !          !  +-----------------31        
##   !         !  !  !          !  !  
##   !         !  !  !          +-29     +----------------30        
##   !         !  !  !             !  +-28  
##   !         !  !  !             !  !  +------------------20        
##   !         !  !  !             +-18  
##   !         !  !  !                !  +----------------12        
##   !         !  !  !                +-10  
##   !         +--8  !                   +----------------10        
##   !            !  !  
##   !            !  ! +------------------13        
##   !            !  +-5 
##   !            !    ! +------------------8         
##   !            !    +-6 
##   !            !      !           +-----------40        
##   !            !      +----------38  
##   !            !                  !  +----------45        
##   !            !                  +-43  
##   !            !                     +----------7         
##   !            ! 
##   !            !                          +37        
##   !            +-------------------------35  
##   !                                       +23        
##   ! 
##   +----------------1         
## 
## 
## remember: this is an unrooted tree!
## 
## Sum of squares =     0.69815
## 
## Average percent standard deviation =     1.87872
## 
## Between        And            Length
## -------        ---            ------
##    1          2                 0.54750
##    1             3              0.03367
##    3            31              0.02399
##   31            32              0.32605
##   32          34                0.38631
##   32          33                0.37573
##   31            26              0.11750
##   26          28                0.52312
##   26             2              0.59233
##    2          4                 0.10576
##    2          3                 0.10870
##    3             4              0.00864
##    4          5                 0.57619
##    4            30              0.02456
##   30            33              0.01375
##   33          35                0.59062
##   33          32                0.58440
##   30            14              0.06373
##   14            42              0.09080
##   42          44                0.58146
##   42            17              0.01512
##   17            22              0.10432
##   22          24                0.49663
##   22          19                0.50322
##   17            16              0.02736
##   16          18                0.55569
##   16            23              0.08971
##   23          25                0.49004
##   23          16                0.46890
##   14             8              0.08143
##    8            21              0.01513
##   21            13              0.01382
##   13            25              0.02493
##   25          27                0.63348
##   25            24              0.00873
##   24          26                0.64713
##   24            37              0.01018
##   37          39                0.79762
##   37            20              0.05733
##   20          22                0.63561
##   20             7              0.01076
##    7          9                 0.60294
##    7            15              0.02045
##   15            27              0.01919
##   27          29                0.62316
##   27          17                0.59833
##   15            12              0.07844
##   12            19              0.01988
##   19          21                0.51958
##   19            41              0.06202
##   41          43                0.41952
##   41          14                0.44222
##   12          6                 0.61125
##   13            11              0.00285
##   11            36              0.03726
##   36          38                0.63767
##   36            40              0.01303
##   40          42                0.61162
##   40          11                0.61742
##   11            39              0.00702
##   39          41                0.63588
##   39             9              0.01201
##    9            34              0.15402
##   34          36                0.52912
##   34          15                0.52914
##    9            29              0.03010
##   29          31                0.60447
##   29            18              0.04380
##   18            28              0.00926
##   28          30                0.56395
##   28          20                0.61623
##   18            10              0.03389
##   10          12                0.54179
##   10          10                0.54035
##   21             5              0.00497
##    5          13                0.63850
##    5             6              0.01245
##    6          8                 0.61868
##    6            38              0.39776
##   38          40                0.39945
##   38            43              0.02264
##   43          45                0.35008
##   43          7                 0.37182
##    8            35              0.91085
##   35          37                0.00025
##   35          23                0.00021
##    1          1                 0.54889
## 
## 
## [1] "Renaming edges"
## [1] "Saving newick data"
## [1] "1291711501383539414412925410"
## [1] "Done:  sliding.rsquared.over.models.output/r.0.73.corecut.1.newick"
## [1] "Doing:  sliding.rsquared.over.models.output/r.0.73.corecut.2.newick"
## [1] "Clustering list"
## [1] "Merging clustering list"
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## [1] "Extracting cores"
## [1] "Extracting MLEs"
```

```
## Calculating MLE for genome 8
## Calculating MLE for genome 38435
## Calculating MLE for genome 68603
## Calculating MLE for genome 74671
## Calculating MLE for genome 81279
## Calculating MLE for genome 139217
## Calculating MLE for genome 247626
## Calculating MLE for genome 953029
## Calculating MLE for genome 1087433
## Calculating MLE for genome 1798918
## Calculating MLE for genome 1864829
## Calculating MLE for genome 1898956
## Calculating MLE for genome 1961502
## Calculating MLE for genome 2233078
## Calculating MLE for genome 4328604
## Calculating MLE for genome 4350340
## Calculating MLE for genome 4381525
## Calculating MLE for genome 6713786
## Calculating MLE for genome 6742235
## Calculating MLE for genome 6753817
## Calculating MLE for genome 6935691
## Calculating MLE for genome 7636370
## Calculating MLE for genome 8069566
## Calculating MLE for genome 8268079
## Calculating MLE for genome 8285225
## Calculating MLE for genome 8298835
## Calculating MLE for genome 8766415
## Calculating MLE for genome 9094481
## Calculating MLE for genome 9105946
## Calculating MLE for genome 10115503
## Calculating MLE for genome 10613611
## Calculating MLE for genome 11252396
## Calculating MLE for genome 11281406
## Calculating MLE for genome 11283907
## Calculating MLE for genome 11287593
## Calculating MLE for genome 11339100
## Calculating MLE for genome 11358995
## Calculating MLE for genome 11557810
## Calculating MLE for genome 11664862
## Calculating MLE for genome 11772215
## Calculating MLE for genome 12124931
## Calculating MLE for genome 12294189
## Calculating MLE for genome 12430564
## Calculating MLE for genome 16244488
## Calculating MLE for genome 16264321
```

```
## [1] "Extracting gprimes"
## [1] "Calculating distance"
## [1] "Running fitch"
```

```
## arguments 'show.output.on.console', 'minimized' and 'invisible' are for Windows only
```

```
## 
##   45 Populations
## 
## Fitch-Margoliash method version 3.696
## 
##                   __ __             2
##                   \  \   (Obs - Exp)
## Sum of squares =  /_ /_  ------------
##                                 2
##                    i  j      Obs
## 
## Negative branch lengths not allowed
## 
## 
##   +----------------2         
##   ! 
##   !              +-----------34        
##   !    +--------32  
##   !    !         +-----------33        
##   ! +-31  
##   ! !  !   +---------------28        
##   ! !  +--26  
##   ! !      !                 +--4         
##   ! !      +-----------------2 
##   ! !                        +---3         
##   ! ! 
##   ! ! +----------------5         
##   1-3 ! 
##   ! ! !     +-----------------35        
##   ! ! !  +-33  
##   ! ! !  !  +-----------------32        
##   ! ! !  !  
##   ! ! !  !     +-----------------44        
##   ! ! !  !     !  
##   ! ! !  !  +-42     +--------------24        
##   ! +-4  !  !  !  +-22  
##   !   !  !  !  !  !  +---------------19        
##   !   !  !  !  +-17  
##   !   !  !  !     !  +----------------18        
##   !   !  !  !     +-16  
##   !   !  !  !        !  +--------------25        
##   !   !  !  !        +-23  
##   !   !  !  !           +-------------16        
##   !   !  !  !  
##   !   +-30  !           +------------------27        
##   !      !  !           !  
##   !      !  !        +-25  +-------------------26        
##   !      !  !        !  !  !  
##   !      !  !        !  +-24  +-----------------------39        
##   !      !  !        !     !  !  
##   !      !  !        !     +-37  +------------------22        
##   !      !  !        !        !  !  
##   !      !  !        !        +-20 +-----------------9         
##   !      !  !        !           ! ! 
##   !      !  !        !           ! !     +------------------29        
##   !      !  !        !           +-7  +-27  
##   !      !  !        !             !  !  +-----------------17        
##   !      !  !        !             !  !  
##   !      !  !     +-13             +-15     +---------------21        
##   !      !  !     !  !                !  +-19  
##   !      +-14     !  !                !  !  !  +------------43        
##   !         !     !  !                +-12  +-41  
##   !         !     !  !                   !     +-------------14        
##   !         !     !  !                   !  
##   !         !     !  !                   +------------------6         
##   !         !     !  !  
##   !         !     !  !     +------------------38        
##   !         !     !  !  +-36  
##   !         !     !  !  !  !  +-----------------42        
##   !         !     !  !  !  +-40  
##   !         !     !  !  !     +------------------11        
##   !         !     !  +-11  
##   !         !     !     !  +------------------41        
##   !         !     !     !  !  
##   !         !  +-21     +-39      +---------------36        
##   !         !  !  !        ! +---34  
##   !         !  !  !        ! !    +---------------15        
##   !         !  !  !        +-9 
##   !         !  !  !          !  +-----------------31        
##   !         !  !  !          !  !  
##   !         !  !  !          +-29     +----------------30        
##   !         !  !  !             !  +-28  
##   !         !  !  !             !  !  +------------------20        
##   !         !  !  !             +-18  
##   !         !  !  !                !  +----------------12        
##   !         !  !  !                +-10  
##   !         +--8  !                   +----------------10        
##   !            !  !  
##   !            !  ! +------------------13        
##   !            !  +-5 
##   !            !    ! +------------------8         
##   !            !    +-6 
##   !            !      !           +-----------40        
##   !            !      +----------38  
##   !            !                  !  +----------45        
##   !            !                  +-43  
##   !            !                     +----------7         
##   !            ! 
##   !            !                          +37        
##   !            +-------------------------35  
##   !                                       +23        
##   ! 
##   +----------------1         
## 
## 
## remember: this is an unrooted tree!
## 
## Sum of squares =     0.69815
## 
## Average percent standard deviation =     1.87872
## 
## Between        And            Length
## -------        ---            ------
##    1          2                 0.54750
##    1             3              0.03367
##    3            31              0.02399
##   31            32              0.32605
##   32          34                0.38631
##   32          33                0.37573
##   31            26              0.11750
##   26          28                0.52312
##   26             2              0.59233
##    2          4                 0.10576
##    2          3                 0.10870
##    3             4              0.00864
##    4          5                 0.57619
##    4            30              0.02456
##   30            33              0.01375
##   33          35                0.59062
##   33          32                0.58440
##   30            14              0.06373
##   14            42              0.09080
##   42          44                0.58146
##   42            17              0.01512
##   17            22              0.10432
##   22          24                0.49663
##   22          19                0.50322
##   17            16              0.02736
##   16          18                0.55569
##   16            23              0.08971
##   23          25                0.49004
##   23          16                0.46890
##   14             8              0.08143
##    8            21              0.01513
##   21            13              0.01382
##   13            25              0.02493
##   25          27                0.63348
##   25            24              0.00873
##   24          26                0.64713
##   24            37              0.01018
##   37          39                0.79762
##   37            20              0.05733
##   20          22                0.63561
##   20             7              0.01076
##    7          9                 0.60294
##    7            15              0.02045
##   15            27              0.01919
##   27          29                0.62316
##   27          17                0.59833
##   15            12              0.07844
##   12            19              0.01988
##   19          21                0.51958
##   19            41              0.06202
##   41          43                0.41952
##   41          14                0.44222
##   12          6                 0.61125
##   13            11              0.00285
##   11            36              0.03726
##   36          38                0.63767
##   36            40              0.01303
##   40          42                0.61162
##   40          11                0.61742
##   11            39              0.00702
##   39          41                0.63588
##   39             9              0.01201
##    9            34              0.15402
##   34          36                0.52912
##   34          15                0.52914
##    9            29              0.03010
##   29          31                0.60447
##   29            18              0.04380
##   18            28              0.00926
##   28          30                0.56395
##   28          20                0.61623
##   18            10              0.03389
##   10          12                0.54179
##   10          10                0.54035
##   21             5              0.00497
##    5          13                0.63850
##    5             6              0.01245
##    6          8                 0.61868
##    6            38              0.39776
##   38          40                0.39945
##   38            43              0.02264
##   43          45                0.35008
##   43          7                 0.37182
##    8            35              0.91085
##   35          37                0.00025
##   35          23                0.00021
##    1          1                 0.54889
## 
## 
## [1] "Renaming edges"
## [1] "Saving newick data"
## [1] "1291711501383603214412927040"
## [1] "Done:  sliding.rsquared.over.models.output/r.0.73.corecut.2.newick"
## [1] "Doing:  sliding.rsquared.over.models.output/r.0.73.corecut.3.newick"
## [1] "Clustering list"
## [1] "Merging clustering list"
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```
## [1] "Extracting cores"
## [1] "Extracting MLEs"
```

```
## Calculating MLE for genome 8
## Calculating MLE for genome 38435
## Calculating MLE for genome 68603
## Calculating MLE for genome 74671
## Calculating MLE for genome 81279
## Calculating MLE for genome 139217
## Calculating MLE for genome 247626
## Calculating MLE for genome 953029
## Calculating MLE for genome 1087433
## Calculating MLE for genome 1798918
## Calculating MLE for genome 1864829
## Calculating MLE for genome 1898956
## Calculating MLE for genome 1961502
## Calculating MLE for genome 2233078
## Calculating MLE for genome 4328604
## Calculating MLE for genome 4350340
## Calculating MLE for genome 4381525
## Calculating MLE for genome 6713786
## Calculating MLE for genome 6742235
## Calculating MLE for genome 6753817
## Calculating MLE for genome 6935691
## Calculating MLE for genome 7636370
## Calculating MLE for genome 8069566
## Calculating MLE for genome 8268079
## Calculating MLE for genome 8285225
## Calculating MLE for genome 8298835
## Calculating MLE for genome 8766415
## Calculating MLE for genome 9094481
## Calculating MLE for genome 9105946
## Calculating MLE for genome 10115503
## Calculating MLE for genome 10613611
## Calculating MLE for genome 11252396
## Calculating MLE for genome 11281406
## Calculating MLE for genome 11283907
## Calculating MLE for genome 11287593
## Calculating MLE for genome 11339100
## Calculating MLE for genome 11358995
## Calculating MLE for genome 11557810
## Calculating MLE for genome 11664862
## Calculating MLE for genome 11772215
## Calculating MLE for genome 12124931
## Calculating MLE for genome 12294189
## Calculating MLE for genome 12430564
## Calculating MLE for genome 16244488
## Calculating MLE for genome 16264321
```

```
## [1] "Extracting gprimes"
## [1] "Calculating distance"
## [1] "Running fitch"
```

```
## arguments 'show.output.on.console', 'minimized' and 'invisible' are for Windows only
```

```
## 
##   45 Populations
## 
## Fitch-Margoliash method version 3.696
## 
##                   __ __             2
##                   \  \   (Obs - Exp)
## Sum of squares =  /_ /_  ------------
##                                 2
##                    i  j      Obs
## 
## Negative branch lengths not allowed
## 
## 
##   +----------------2         
##   ! 
##   !              +-----------34        
##   !    +--------32  
##   !    !         +-----------33        
##   ! +-31  
##   ! !  !   +---------------28        
##   ! !  +--26  
##   ! !      !                 +--4         
##   ! !      +-----------------2 
##   ! !                        +---3         
##   ! ! 
##   ! ! +----------------5         
##   1-3 ! 
##   ! ! !     +-----------------35        
##   ! ! !  +-33  
##   ! ! !  !  +-----------------32        
##   ! ! !  !  
##   ! ! !  !     +-----------------44        
##   ! ! !  !     !  
##   ! ! !  !  +-42     +--------------24        
##   ! +-4  !  !  !  +-22  
##   !   !  !  !  !  !  +---------------19        
##   !   !  !  !  +-17  
##   !   !  !  !     !  +----------------18        
##   !   !  !  !     +-16  
##   !   !  !  !        !  +--------------25        
##   !   !  !  !        +-23  
##   !   !  !  !           +-------------16        
##   !   !  !  !  
##   !   +-30  !           +------------------27        
##   !      !  !           !  
##   !      !  !        +-25  +-------------------26        
##   !      !  !        !  !  !  
##   !      !  !        !  +-24  +-----------------------39        
##   !      !  !        !     !  !  
##   !      !  !        !     +-37  +------------------22        
##   !      !  !        !        !  !  
##   !      !  !        !        +-20 +-----------------9         
##   !      !  !        !           ! ! 
##   !      !  !        !           ! !     +------------------29        
##   !      !  !        !           +-7  +-27  
##   !      !  !        !             !  !  +-----------------17        
##   !      !  !        !             !  !  
##   !      !  !     +-13             +-15     +---------------21        
##   !      !  !     !  !                !  +-19  
##   !      +-14     !  !                !  !  !  +------------43        
##   !         !     !  !                +-12  +-41  
##   !         !     !  !                   !     +-------------14        
##   !         !     !  !                   !  
##   !         !     !  !                   +------------------6         
##   !         !     !  !  
##   !         !     !  !     +------------------38        
##   !         !     !  !  +-36  
##   !         !     !  !  !  !  +-----------------42        
##   !         !     !  !  !  +-40  
##   !         !     !  !  !     +------------------11        
##   !         !     !  +-11  
##   !         !     !     !  +------------------41        
##   !         !     !     !  !  
##   !         !  +-21     +-39      +---------------36        
##   !         !  !  !        ! +---34  
##   !         !  !  !        ! !    +---------------15        
##   !         !  !  !        +-9 
##   !         !  !  !          !  +-----------------31        
##   !         !  !  !          !  !  
##   !         !  !  !          +-29     +----------------30        
##   !         !  !  !             !  +-28  
##   !         !  !  !             !  !  +------------------20        
##   !         !  !  !             +-18  
##   !         !  !  !                !  +----------------12        
##   !         !  !  !                +-10  
##   !         +--8  !                   +----------------10        
##   !            !  !  
##   !            !  ! +------------------13        
##   !            !  +-5 
##   !            !    ! +------------------8         
##   !            !    +-6 
##   !            !      !           +-----------40        
##   !            !      +----------38  
##   !            !                  !  +----------45        
##   !            !                  +-43  
##   !            !                     +----------7         
##   !            ! 
##   !            !                          +37        
##   !            +-------------------------35  
##   !                                       +23        
##   ! 
##   +----------------1         
## 
## 
## remember: this is an unrooted tree!
## 
## Sum of squares =     0.69815
## 
## Average percent standard deviation =     1.87872
## 
## Between        And            Length
## -------        ---            ------
##    1          2                 0.54750
##    1             3              0.03367
##    3            31              0.02399
##   31            32              0.32605
##   32          34                0.38631
##   32          33                0.37573
##   31            26              0.11750
##   26          28                0.52312
##   26             2              0.59233
##    2          4                 0.10576
##    2          3                 0.10870
##    3             4              0.00864
##    4          5                 0.57619
##    4            30              0.02456
##   30            33              0.01375
##   33          35                0.59062
##   33          32                0.58440
##   30            14              0.06373
##   14            42              0.09080
##   42          44                0.58146
##   42            17              0.01512
##   17            22              0.10432
##   22          24                0.49663
##   22          19                0.50322
##   17            16              0.02736
##   16          18                0.55569
##   16            23              0.08971
##   23          25                0.49004
##   23          16                0.46890
##   14             8              0.08143
##    8            21              0.01513
##   21            13              0.01382
##   13            25              0.02493
##   25          27                0.63348
##   25            24              0.00873
##   24          26                0.64713
##   24            37              0.01018
##   37          39                0.79762
##   37            20              0.05733
##   20          22                0.63561
##   20             7              0.01076
##    7          9                 0.60294
##    7            15              0.02045
##   15            27              0.01919
##   27          29                0.62316
##   27          17                0.59833
##   15            12              0.07844
##   12            19              0.01988
##   19          21                0.51958
##   19            41              0.06202
##   41          43                0.41952
##   41          14                0.44222
##   12          6                 0.61125
##   13            11              0.00285
##   11            36              0.03726
##   36          38                0.63767
##   36            40              0.01303
##   40          42                0.61162
##   40          11                0.61742
##   11            39              0.00702
##   39          41                0.63588
##   39             9              0.01201
##    9            34              0.15402
##   34          36                0.52912
##   34          15                0.52914
##    9            29              0.03010
##   29          31                0.60447
##   29            18              0.04380
##   18            28              0.00926
##   28          30                0.56395
##   28          20                0.61623
##   18            10              0.03389
##   10          12                0.54179
##   10          10                0.54035
##   21             5              0.00497
##    5          13                0.63850
##    5             6              0.01245
##    6          8                 0.61868
##    6            38              0.39776
##   38          40                0.39945
##   38            43              0.02264
##   43          45                0.35008
##   43          7                 0.37182
##    8            35              0.91085
##   35          37                0.00025
##   35          23                0.00021
##    1          1                 0.54889
## 
## 
## [1] "Renaming edges"
## [1] "Saving newick data"
## [1] "1291711501383651514412928640"
## [1] "Done:  sliding.rsquared.over.models.output/r.0.73.corecut.3.newick"
```
