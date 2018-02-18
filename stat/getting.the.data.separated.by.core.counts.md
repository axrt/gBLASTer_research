# Getting the data separated by core counts

First to get the libraries.


```r
if(!require("gbra")){
  install.packages("gbra")
  library("gbra")
}
```

```
## Loading required package: gbra
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
## The following objects are masked from 'package:stats':
## 
##     filter, lag
## 
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
## 
## Loading required package: data.table
## 
## Attaching package: 'data.table'
## 
## The following objects are masked from 'package:dplyr':
## 
##     between, last
## 
## Loading required package: tidyr
## Loading required package: Rcpp
## Loading required package: boot
## Loading required package: stringr
```

```r
if(!require("reshape2")){
  install.packages("reshape2")
  library("reshape2")
}
```

```
## Loading required package: reshape2
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
#Loading data

Load the data (in case not there yet)


```r
bh.data <- read.table(file = "bh.data.minhit_4.txt")
legend<-load.legend(legend.file = "/home/alext/Developer/gbra/inst/extradata/legend.csv",sep = "\t")
```

#Data processing

Processing the data to split it into three groups:

* the one that only contains core genes
* the one that has outsiders just once
* the rest of the genes, that can be outsiders much more frequently


```r
dir<-"data" #data folder
factor.core<-c("core","outsider","nosep") #generate three factors to mark states
if(!file.exists("master.table.core.txt")){
  #generate a matrix big enough to hold the data (essentially the same size as the bh.data), initially filled in with "outsider"
  master.table.core<-as.data.frame(matrix(nrow=nrow(bh.data),ncol=nrow(legend),factor.core[3]),stringsAsFactors = FALSE)
  #copy the row names from bh.data
  row.names(master.table.core)<-row.names(bh.data)
  #copy colnames fromt he legend
  colnames(master.table.core)<-legend$id_genomes
  #now, for each row in the initial bh.table go over the models one by one, load them into ram and see where the row gene is placed 
  #relative to each genome
  for(i in 1:(nrow(legend)-1)){
    #j must be i+1 because the diagonas are always "nonsep" as you can't really separate the organism with itself
    for(j in (i+1):nrow(legend)){
      #this if is actually redundant, need to remove
      if(i!=j){
        #load the corresponding model
        load(file = paste(dir,paste(paste(legend$id_genomes[i],legend$id_genomes[j],sep="_"),".rda",sep = ""),sep = "/"))
        #see if the model R squared score if over 0.75, which makes sence in terms of true positive ratio (needs link)
        if(as.numeric(l.m.g$Rsq)>=0.75){
          #for the pair of genomes - restrict the bh.table to only those selected
          data.table<-select.genomes(df = bh.data, g.ids = c(legend$id_genomes[i],legend$id_genomes[j]))
          #attempt to predict the probabilities for each pair with a precalculated logistic regression
          glm.probs<-predict(l.m.g$fit, data.table, type = "response")
          #actually split into two groups with a 0.5 probability cutoff
          glm.probs.groups<-sapply(glm.probs,function(i){if(i>=0.5){return(FALSE)}else{return(TRUE)}})
          #append another column, which can be either true of false to indicate the clusters 
          data.table$cluster<-grepl(pattern = paste("\\b",legend$id_genomes[i],"X",sep = "",collapse = ""),x = rownames(data.table),perl = TRUE)
          #see if the clustering assignment matches the true clusters, mark as core if yes, otherwise - outsider
          sep.factor.v<-ifelse(test = data.table$cluster==glm.probs.groups, yes = factor.core[1], no = factor.core[2])
          #assign common name identificators like "X1,2..." in accordance with the legend
          sep.factor.v.i<-sep.factor.v[grepl(pattern = paste("\\b",legend$id_genomes[i],"X",sep = "",collapse = ""),x = names(sep.factor.v),perl = TRUE)]
          sep.factor.v.j<-sep.factor.v[grepl(pattern = paste("\\b",legend$id_genomes[j],"X",sep = "",collapse = ""),x = names(sep.factor.v),perl = TRUE)]
          core.col.j<-which(legend$id_genomes==legend$id_genomes[j])
          core.col.i<-which(legend$id_genomes==legend$id_genomes[i])
          #mark i and j cells (because the matrix here is not triangle and needs symmetry)
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
  #finally save out the resulting table
  write.table(x=master.table.core, file="master.table.core.txt", sep="\t", quote = FALSE)
}
```

#Gene states evaluation

For each genomes calculate the totals of states for the genes


```r
#read the data
master.table.core<-read.table(file="master.table.core.txt")
#create a parallel cluster with the cluster log output
file.remove("cluster.txt")
```

```
## [1] TRUE
```

```r
cluster <<- makeCluster(6,outfile = "cluster.txt")
registerDoParallel(cluster)
#create a diriving summary of genes
master.table.core.per.gene.summary<-as.data.frame(t(parApply(cl=cluster, X = master.table.core, MARGIN = 1, FUN = function(i){
  sum.core<-length(i[i=="core"])
  sum.outsider<-length(i[i=="outsider"])
  sum.nosep<-length(i[i=="nosep"])
  return(c(core=sum.core, outsider=sum.outsider, nonsep=sum.nosep))
})))
#append the names because the following filter will erase the rownames
master.table.core.per.gene.summary$gene.name<-row.names(master.table.core.per.gene.summary)
#do not forget to cluse the cluster
stopCluster(cluster)
#separate the "always core table"
master.table.core.always.core<- master.table.core.per.gene.summary %>% filter(outsider==0)
head(master.table.core.always.core)
```

```
##   core outsider nonsep gene.name
## 1   15        0     30       1X8
## 2   15        0     30      1X17
## 3   15        0     30      1X59
## 4   15        0     30      1X70
## 5   15        0     30      1X82
## 6   15        0     30      1X91
```

```r
#separate the "outsider once"
master.table.core.outsider.once<- master.table.core.per.gene.summary %>% filter(outsider==1)
head(master.table.core.outsider.once)
```

```
##   core outsider nonsep gene.name
## 1   14        1     30      1X54
## 2   14        1     30      1X95
## 3   14        1     30     1X156
## 4   14        1     30     1X679
## 5   14        1     30     1X699
## 6   14        1     30     1X741
```

```r
#separate the frequent outsider
master.table.core.outsider.frequently<- master.table.core.per.gene.summary %>% filter(outsider>1)
head(master.table.core.outsider.frequently)
```

```
##   core outsider nonsep gene.name
## 1   12        3     30      1X10
## 2    6        9     30      1X37
## 3   10        5     30      1X39
## 4   13        2     30     1X101
## 5   12        3     30     1X164
## 6   12        3     30     1X453
```

```r
#write out the data
write.table(master.table.core.always.core, file="master.table.core.always.core.txt", sep="\t", col.names=TRUE, quote = FALSE)
write.table(master.table.core.outsider.once, file="master.table.core.outsider.once.txt", sep="\t", col.names=TRUE, quote = FALSE)
write.table(master.table.core.outsider.frequently, file="master.table.core.outsider.frequently.txt", sep="\t", col.names=TRUE, quote = FALSE)
#
```
