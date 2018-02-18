# Bootstrapped Models
Alexander Tuzhikov  
September 29, 2016  

The idea behind this iteration is: for each pair of organisms we will create 1000 models. From each model we will take the betta coefficients and mean them. Create a resulting model, which will have the coefficients, which are weighted by the number of votes of all models, created from different substes of the data.

```r
#compile gbra to ensure it's up to date
library(devtools)
install_github("axrt/gbra")#actually at this moment there is a lot of functionality that needs to be updated
library(gbra)
source("R/helper.R")
sixtyseven.genomes<-"67.genomes"
model.df.file<- "model.df.txt"
bh.data.raw.file<- "bh.data.raw.rda"
bh.data.normal.file<- "bh.data.normal.rda"
library(doMC)
library(doParallel)
library(foreach)
library(ggplot2)
library(data.table)
Rsq.cutoff<- 0.73
logfile<- "log.txt"
library(plyr)
library(dplyr)
install_github("vqv/ggbiplot")
library(boot)
library(pander)
library(seplyr)
```

We create a legend, that will only contain the genomes that we want to include (the others are not good enough for this or that reason)

```r
legend<- read.table(file=sixtyseven.genomes.dir("legend.txt"), sep="\t", header=TRUE, stringsAsFactors = FALSE) %>% 
  dplyr::filter(id_genomes!=23) %>% 
  dplyr::filter(id_genomes!=40) %>% 
  dplyr::filter(id_genomes!=47) %>%
  dplyr::filter(id_genomes!=57) %>%
  dplyr::filter(id_genomes!=63) %>%
  dplyr::filter(id_genomes!=64) %>%
  dplyr::select(name, id_genomes) %>%
  dplyr::arrange(id_genomes)
```

We load the bidirectional hits form the previous iterations of analaysis.

```r
load(sixtyseven.genomes.dir(bh.data.raw.file))
```

Convert the bh data to a normalized form, genome ids and minimal hits restricted (4 in this case)

```r
bh.data%>% 
    as.data.frame %>%
    restrict.minimal.hits(minhit = 4) %>%
    normalize.scores %>%
    attach.genomeid.header %>%
    sign.bh.table -> bh.data.norm
rm(bh.data) #saves ram quite a bit
bh.data.norm %>% head %>% pander(caption="Best Hit DataFrame")
```


-------------------------------------------------------------------------------
  &nbsp;    X1     X10      X11      X12      X13      X14      X15      X16   
---------- ---- --------- -------- -------- -------- -------- -------- --------
 **1X8**    1    0.09708     0      0.117    0.1585     0        0      0.1388 

 **1X10**   1       0        0      0.1247     0        0        0        0    

 **1X17**   1    0.2406      0      0.2153   0.2231     0        0        0    

 **1X37**   1    0.4631      0      0.4824   0.4378   0.4515     0        0    

 **1X39**   1       0        0        0      0.2079     0        0        0    

 **1X54**   1    0.5643    0.5771   0.5455   0.5721   0.5354   0.4596   0.5231 
-------------------------------------------------------------------------------

Table: Best Hit DataFrame (continued below)

 
------------------------------------------------------------------------
  &nbsp;     X17       X18     X19     X2       X20      X21      X22   
---------- -------- --------- ----- --------- -------- -------- --------
 **1X8**      0      0.04753    0    0.08968   0.063      0        0    

 **1X10**     0      0.1484     0     0.155      0        0        0    

 **1X17**   0.1214      0       0       0      0.1624     0        0    

 **1X37**   0.1299      0       0       0      0.1531   0.4141   0.1498 

 **1X39**     0         0       0       0        0      0.2116     0    

 **1X54**   0.4303   0.4503     0    0.5291    0.4121   0.436    0.578  
------------------------------------------------------------------------

Table: Table continues below

 
-----------------------------------------------------------------------------
  &nbsp;      X24       X25      X26      X27       X28       X3       X30   
---------- --------- --------- -------- -------- --------- --------- --------
 **1X8**    0.05268   0.05515     0      0.1047   0.09976   0.06434     0    

 **1X10**   0.1366       0        0        0         0         0        0    

 **1X17**      0         0        0      0.1594   0.2742       0        0    

 **1X37**      0         0        0      0.522    0.1151       0        0    

 **1X39**      0         0        0        0         0         0        0    

 **1X54**      0      0.4769    0.5947   0.5745   0.6643    0.6222    0.4321 
-----------------------------------------------------------------------------

Table: Table continues below

 
--------------------------------------------------------------------------
  &nbsp;     X31      X32      X33      X34      X35       X36      X37   
---------- -------- -------- -------- -------- -------- --------- --------
 **1X8**    0.1269   0.1069   0.1004     0        0      0.06995     0    

 **1X10**     0        0      0.1247     0        0         0        0    

 **1X17**   0.1205   0.2476     0        0        0         0        0    

 **1X37**   0.3012   0.5424     0        0        0         0        0    

 **1X39**     0        0        0        0        0         0        0    

 **1X54**   0.533    0.4551   0.5536   0.4488   0.2733   0.5058    0.3969 
--------------------------------------------------------------------------

Table: Table continues below

 
----------------------------------------------------------------------------
  &nbsp;     X38      X39      X4        X41       X42       X43      X44   
---------- -------- ------- --------- --------- --------- --------- --------
 **1X8**    0.2901     0     0.05672   0.08945   0.05044   0.07353     0    

 **1X10**     0        0        0      0.1669       0      0.1276      0    

 **1X17**     0        0        0         0      0.1672     0.231      0    

 **1X37**     0        0        0         0      0.2285    0.3469    0.5165 

 **1X39**   0.2033     0        0         0         0         0        0    

 **1X54**   0.6228   0.567   0.6329    0.4193     0.442    0.4688    0.459  
----------------------------------------------------------------------------

Table: Table continues below

 
---------------------------------------------------------------------------
  &nbsp;     X45      X46       X48       X49       X5      X50      X51   
---------- -------- -------- --------- --------- -------- -------- --------
 **1X8**    0.1007     0      0.07017   0.08788   0.1004   0.1592     0    

 **1X10**     0      0.1657   0.2244    0.1651      0      0.1306   0.1318 

 **1X17**   0.2956     0         0         0      0.2498     0        0    

 **1X37**     0        0      0.5011    0.5055      0      0.1239     0    

 **1X39**     0        0         0         0        0        0        0    

 **1X54**   0.5076   0.4264   0.5679    0.5386    0.6258   0.4515     0    
---------------------------------------------------------------------------

Table: Table continues below

 
--------------------------------------------------------------------------------
  &nbsp;      X52      X53      X54      X55      X56     X58    X59       X6   
---------- --------- -------- -------- -------- -------- ----- -------- --------
 **1X8**    0.07937     0        0        0        0       0      0        0    

 **1X10**   0.1948      0        0        0        0       0      0        0    

 **1X17**      0        0        0        0      0.1157    0      0        0    

 **1X37**      0      0.1944   0.4708     0        0       0      0      0.4477 

 **1X39**   0.4787      0      0.1987     0        0       0      0        0    

 **1X54**   0.5189    0.3954   0.4933   0.4509   0.4414    0    0.4124     0    
--------------------------------------------------------------------------------

Table: Table continues below

 
-------------------------------------------------------------------------
  &nbsp;     X60       X62      X65      X66      X67      X68      X69  
---------- -------- --------- -------- -------- -------- -------- -------
 **1X8**      0      0.08429     0        0      0.2874     0        0   

 **1X10**     0         0        0      0.1336     0        0        0   

 **1X17**     0         0      0.1581   0.1554   0.1716     0        0   

 **1X37**     0         0      0.5187   0.5154   0.1239   0.4124     0   

 **1X39**     0         0        0        0      0.2375     0        0   

 **1X54**   0.3154   0.6368    0.4619   0.5428   0.5837   0.4375   0.444 
-------------------------------------------------------------------------

Table: Table continues below

 
-------------------------------------
  &nbsp;      X7       X8       X9   
---------- -------- -------- --------
 **1X8**      0      0.1455     0    

 **1X10**   0.1645   0.1259     0    

 **1X17**     0      0.1118     0    

 **1X37**     0      0.1399   0.2197 

 **1X39**     0        0        0    

 **1X54**   0.4187   0.6117   0.4187 
-------------------------------------

Now, create the model driver data frame, which will contain model pairs (see below).

```r
model.driver<- create.model.driver(legend)
model.driver %>% head %>% pander(caption="Model Driver DataFrame")
```


-------------------------
 X1   X2       file      
---- ---- ---------------
 1    2    1_2_model.rda 

 1    3    1_3_model.rda 

 1    4    1_4_model.rda 

 1    5    1_5_model.rda 

 1    6    1_6_model.rda 

 1    7    1_7_model.rda 
-------------------------

Table: Model Driver DataFrame

Create a log file for the multitheaded processing below, remove the file if it exists.

```r
model.training.log.file<- "model.training.log.txt" %>% sixtyseven.genomes.quadruple.model
if(model.training.log.file %>% file.exists){
        file.remove(model.training.log.file)
}
```

`logreg.model.bootstrap` is going to be a funciton that will return betta coefficients for models that we are going to get in the bootstrap process.

```r
logreg.model.bootstrap<-function(data, indices, formula){
  working.data<-data[indices,]
  glm.fit<-glm(formula = formula, data=working.data, family="binomial")
  return(glm.fit$coefficients)
}
```

Train models, skip if the model has already been trained.

```r
models<- foreach(i=1:nrow(model.driver)) %do%{
  
  org.pair<- c(model.driver$X1[i], model.driver$X2[i])
  out.file.stepped<- sixtyseven.genomes.quadruple.model.bootstrapped(paste0(org.pair[1], "_", org.pair[2],".stepped.rds"))
  out.file.bootstrapped<- sixtyseven.genomes.quadruple.model.bootstrapped(paste0(org.pair[1], "_", org.pair[2],".bootstrapped.rds"))
  if(file.exists(out.file.bootstrapped)){
    return(out.file.bootstrapped)
  }
  message(paste("Doing", org.pair[1], "::", org.pair[2]))
  
  coocked.data<- bh.data.norm %>% select.genome.data(org.pair) %>%
      filter.genome.data(org.pair) %>% 
      append.cluster %>% data.frame
  
  f<- get.glm.formula(coocked.data)
  glm.fit<- glm(formula = f, data=coocked.data, family="binomial")
  glm.fit<- step(glm.fit)
  saveRDS(object = glm.fit, 
          file = out.file.stepped,
          compress = "gzip")
  model<-boot(glm.fit$model, logreg.model.bootstrap, R=1000, formula=glm.fit$formula, parallel = "multicore", ncpus = 6)
  model.coeff<- apply(model$t, 2, median)
  names(model.coeff)<- names(model$t0)
  glm.fit$coefficients<- model.coeff
  
  saveRDS(object = glm.fit, 
          file = out.file.bootstrapped,
          compress = "gzip")
  return(out.file.bootstrapped)
}
```

We append paths to the model files to the model driver and save the driver to a file to use later.

```r
registerDoMC(cores=6)
model.driver.file<- "model.driver.bootstrapped.txt" %>% sixtyseven.genomes.quadruple.model
if(!model.driver.file %>% file.exists){
        
    model.driver$rsquared<- foreach(i=1:nrow(model.driver), .combine="c") %dopar%{
    org.pair<- c(model.driver$X1[i], model.driver$X2[i])
    out.file<- sixtyseven.genomes.quadruple.model.bootstrapped(paste0(org.pair[1], "_", org.pair[2],".bootstrapped.rds"))
   model<- readRDS(file=out.file)
  
   return(get.rsq(model))
    }
    write.table(model.driver, file=model.driver.file, sep="\t", quote = FALSE)
}else{
    model.driver<- read.table(model.driver.file, sep="\t", header=TRUE, stringsAsFactors = FALSE)    
}
model.driver %>% head %>% pander(caption="Updated Model Driver DataFrame")
```


------------------------------------
 X1   X2       file        rsquared 
---- ---- --------------- ----------
 1    2    1_2_model.rda   0.05169  

 1    3    1_3_model.rda    0.8749  

 1    4    1_4_model.rda    0.8742  

 1    5    1_5_model.rda   0.08929  

 1    6    1_6_model.rda    0.9173  

 1    7    1_7_model.rda    0.8528  
------------------------------------

Table: Updated Model Driver DataFrame


Now we will go over all models and collect the R-squared values. That will be for *D.simulans*, and we will plot a line at R-squared = 0.73 because this value roughly corresponds to 95% specificity.

```r
rsq.plot<- model.driver %>%
  rename(query.genome=X1, target.genome=X2) %>% 
  filter(query.genome==6|target.genome==6) %>%
  mutate(target.genome=sapply(1:length(.$target.genome), function(x){
  if(.$target.genome[x]==6){
    return(.$query.genome[x])
  }else{
    return(.$target.genome[x])
  }
})) %>%
  select(rsquared, target.genome) %>%
  merge(x=., y=legend, by.x="target.genome", by.y="id_genomes") %>%
  mutate(name=str_trim(toupper(str_replace_all(string=.$name, pattern = "\\.fasta|_|complete_genome|chromosome", replacement=" ")))) %>% 
  arrange(rsquared) %>%
  mutate(name=factor(.$name, levels=.$name)) %>%
  melt(id.vars="name") %>%
  filter(variable %in% c("rsquared")) %>%
  ggplot(data=., mapping=aes(x=name, y=value, group=variable, color=variable)) + geom_point() +
  geom_smooth() +
  theme_bw() + theme(axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none", plot.margin = unit(c(1, 1, 3, 2), "cm"), axis.title.x=element_text(vjust=1)) + 
  geom_hline(yintercept = 0.73, color="red") + 
  xlab("Genome") + ylab("R-squared") + labs(title = "Drosophila simulans models R-squared values")

plot(rsq.plot)
```

![](bootstrapped.models_files/figure-html/rsquared plot for Drosophila-1.png)<!-- -->

```r
pdf(file = "rsq.drosoph.bootstrapped.pdf"%>% sixtyseven.genomes.quadruple.model, width = 15, height=8)
plot(rsq.plot)
dev.off()
```

```
## png 
##   2
```

I have been stealing some code from previous analysis, so here I will have some hacks to have it compatible with the current iteration.
The huuuuge chunck of code below will create a so called master-core-table for our bootstrapped models. It will show for each orf on the analysis for each organism if we could crearly separate the given orf, and, if we could, whether it went to core, or did not and thus fell into the cluster of the opposing genome. We chose wide margin of (0.025, 0.975) of probability which keeps the orf as "not separated" for any orf that had cluster assignment probability fall into this margin. In human language (mostly, if I put myself together to write it better), if we have an orf and a model predicts it's porbability to go into cluster A (from A and B cluster pair) as say 0.6 - it is not good enough and we say, that we do not have enough evidence ot have this orf assigned to a particular cluster. But for a case of probability 0.98, of an orf to be from a certain cluster - we say that yes, this is good enough and this orf is comeing from its origin cluster (core) or from an opposing cluster (outsider).

```r
bh.data<- bh.data.norm
factor.core<-c("core","outsider","nosep") #generate three factors to mark states
if(!"master.table.core.bootstrapped.0025.txt" %>% sixtyseven.genomes.quadruple.model %>% file.exists){
  #generate a matrix big enough to hold the data (essentially the same size as the bh.data), initially filled in with "outsider"
  master.table.core<-as.data.frame(
          matrix(
                  nrow=nrow(bh.data),
                  ncol=nrow(legend),
                  factor.core[3]),
          stringsAsFactors = FALSE)
  #copy the row names from bh.data
  row.names(master.table.core)<-row.names(bh.data)
  #copy colnames from the legend
  colnames(master.table.core)<-legend$id_genomes
  #now, for each row in the initial bh.table go over the models one by one, load them into ram and see where the row gene is placed 
  #relative to each genome
  registerDoMC(cores=6)
  writeLines(c(""), "log.txt")
  master.table.core.list<-foreach(i=1:(nrow(legend)-1))%do%{
    #j must be i+1 because the diagonas are always "nonsep" as you can't really separate the organism with itself
    foreach(j=(i+1):nrow(legend))%dopar%{
        sink("log.txt", append=TRUE)
        #load the corresponding model list
        model.file<- paste(
                paste(
                        legend$id_genomes[i],
                        legend$id_genomes[j],
                        sep="_"),
                ".bootstrapped.rds",
                sep = "") %>% sixtyseven.genomes.quadruple.model.bootstrapped
        l.m.g<- readRDS(file=model.file)
        #see if all models are present in the list (which indicates that the separation can be done definitively in all attempts with a given rsquared)
        if(get.rsq(l.m.g)>=0.73){
                #select corresponding genomes
                data.table<-select.genome.data(bh.data, c(legend$id_genomes[i],legend$id_genomes[j])) %>% filter.genome.data(c(legend$id_genomes[i],legend$id_genomes[j]))
                #attempt to predict the probabilities for each pair with a precalculated logistic regression
                data.table$cluster<-grepl(pattern = paste("\\b",legend$id_genomes[i],"X",sep = "",collapse = ""),x = rownames(data.table),perl = TRUE)
                glm.probs<- predict(l.m.g, data.table, type = "response")
                glm.probs.groups<- sapply(glm.probs,function(l){
                                if(l<0.025){
                                        return(TRUE)
                                }
                                else if(l>0.975){
                                        return(FALSE)
                                }else{
                                      return(NA)
                                }
                                })
                
                sep.factor.v<- sapply(1:length(data.table$cluster), function(i){
                        if(is.na(glm.probs.groups[i])){
                                return(factor.core[3])
                        }
                        if(data.table$cluster[i]==glm.probs.groups[i]){
                                return(factor.core[2])
                        }else{
                                return(factor.core[1])
                        }
                })
                sep.factor.v<- ifelse(test=data.table$cluster==glm.probs.groups, yes=factor.core[2], no=factor.core[1])

                output.table<- cbind(name=rownames(data.table), cluster=data.table$cluster, match=sep.factor.v==factor.core[1]) %>% data.frame
                write.table(x=output.table, file=paste(legend$id_genomes[i], legend$id_genomes[j], sep="_") %>% 
                                    paste0(".txt.bootstrapped.match") %>% sixtyseven.genomes.quadruple.model.match(), 
                            sep="\t", quote = FALSE, row.names = FALSE)
          #assign common name identificators like "X1,2..." in accordance with the legend
          sep.factor.v.i<-sep.factor.v[grepl(pattern = paste("\\b",legend$id_genomes[i],"X",sep = "",collapse = ""),x = names(sep.factor.v),perl = TRUE)]
          sep.factor.v.j<-sep.factor.v[grepl(pattern = paste("\\b",legend$id_genomes[j],"X",sep = "",collapse = ""),x = names(sep.factor.v),perl = TRUE)]
          core.col.j<-which(legend$id_genomes==legend$id_genomes[j])
          core.col.i<-which(legend$id_genomes==legend$id_genomes[i])
          #mark i and j cells (because the matrix here is not triangle and needs symmetry)
          return(list(sep.factor.v.i, sep.factor.v.j, core.col.i, core.col.j))
        }else{
                return(list(i, j))
        }
    }
  }
  for(i in 1:length(master.table.core.list)){
          for(j in 1:length(master.table.core.list[[i]])){
                  if(length(master.table.core.list[[i]][[j]])==4){
                          print(paste("Doing:", i, j, sep=" "))
                          sep.factor.v.i<- master.table.core.list[[i]][[j]][[1]]
                          sep.factor.v.j<- master.table.core.list[[i]][[j]][[2]]
                          core.col.i<- master.table.core.list[[i]][[j]][[3]]
                          core.col.j<- master.table.core.list[[i]][[j]][[4]]
                          master.table.core[names(sep.factor.v.i),core.col.j]<-sep.factor.v.i
                          master.table.core[names(sep.factor.v.j),core.col.i]<-sep.factor.v.j
                  }
          }
  }
  
  #finally save out the resulting table
  write.table(x=master.table.core, file="master.table.core.bootstrapped.0025.txt" %>% sixtyseven.genomes.quadruple.model, sep="\t", quote = FALSE)
}
```


```r
#read the data
master.table.core<-read.table(file="master.table.core.bootstrapped.0025.txt" %>% sixtyseven.genomes.quadruple.model)
master.table.core[is.na(master.table.core)]<- factor.core[3]
master.table.core %>% head %>% pander(caption="Master Table Core DataFrame")
```


--------------------------------------------------------------------------
  &nbsp;     X1      X2      X3      X4      X5      X6      X7      X8   
---------- ------- ------- ------- ------- ------- ------- ------- -------
 **1X8**    nosep   nosep   core    core    nosep   core    nosep   nosep 

 **1X10**   nosep   nosep   core    core    nosep   core    nosep   nosep 

 **1X17**   nosep   nosep   core    core    nosep   core    core    nosep 

 **1X37**   nosep   nosep   core    core    nosep   nosep   core    nosep 

 **1X39**   nosep   nosep   core    core    nosep   core    nosep   nosep 

 **1X54**   nosep   nosep   nosep   nosep   nosep   core    nosep   nosep 
--------------------------------------------------------------------------

Table: Master Table Core DataFrame (continued below)

 
--------------------------------------------------------------------------------
  &nbsp;       X9       X10     X11     X12     X13      X14       X15     X16  
---------- ---------- ------- ------- ------- ------- ---------- ------- -------
 **1X8**     nosep     nosep   nosep   nosep   nosep     core     core    nosep 

 **1X10**    nosep     nosep   nosep   nosep   nosep    nosep     core    nosep 

 **1X17**    nosep     nosep   nosep   nosep   nosep     core     core    nosep 

 **1X37**   outsider   nosep   nosep   nosep   nosep   outsider   nosep   nosep 

 **1X39**    nosep     nosep   nosep   nosep   nosep    nosep     core    nosep 

 **1X54**    nosep     nosep   nosep   nosep   nosep     core     core    nosep 
--------------------------------------------------------------------------------

Table: Table continues below

 
------------------------------------------------------------------------------
  &nbsp;      X17       X18     X19      X20        X21        X22       X24  
---------- ---------- ------- ------- ---------- ---------- ---------- -------
 **1X8**      core     nosep   nosep    nosep       core       core     nosep 

 **1X10**     core     nosep   nosep    nosep       core      nosep     nosep 

 **1X17**     core     nosep   nosep    nosep       core      nosep     nosep 

 **1X37**   outsider   nosep   nosep   outsider   outsider   outsider   nosep 

 **1X39**    nosep     nosep   nosep    nosep      nosep      nosep     nosep 

 **1X54**     core     nosep   nosep     core      nosep      nosep     nosep 
------------------------------------------------------------------------------

Table: Table continues below

 
-----------------------------------------------------------------------------
  &nbsp;     X25     X26     X27     X28      X30       X31     X32     X33  
---------- ------- ------- ------- ------- ---------- ------- ------- -------
 **1X8**    nosep   nosep   nosep   nosep     core     nosep   nosep   nosep 

 **1X10**   nosep   nosep   nosep   nosep    nosep     nosep   nosep   nosep 

 **1X17**   nosep   nosep   nosep   nosep     core     nosep   nosep   nosep 

 **1X37**   nosep   nosep   nosep   nosep   outsider   nosep   nosep   nosep 

 **1X39**   nosep   nosep   nosep   nosep    nosep     nosep   nosep   nosep 

 **1X54**   nosep   nosep   nosep   nosep    nosep     nosep   nosep   nosep 
-----------------------------------------------------------------------------

Table: Table continues below

 
--------------------------------------------------------------------------
  &nbsp;     X34     X35     X36     X37     X38     X39     X41     X42  
---------- ------- ------- ------- ------- ------- ------- ------- -------
 **1X8**    nosep   nosep   nosep   core    nosep   nosep   nosep   nosep 

 **1X10**   nosep   nosep   nosep   nosep   nosep   nosep   nosep   nosep 

 **1X17**   nosep   nosep   nosep   core    nosep   nosep   nosep   nosep 

 **1X37**   nosep   nosep   nosep   nosep   nosep   nosep   nosep   nosep 

 **1X39**   nosep   nosep   nosep   nosep   nosep   nosep   nosep   nosep 

 **1X54**   nosep   nosep   nosep   core    nosep   nosep   nosep   nosep 
--------------------------------------------------------------------------

Table: Table continues below

 
-----------------------------------------------------------------------------
  &nbsp;     X43      X44       X45     X46     X48     X49     X50     X51  
---------- ------- ---------- ------- ------- ------- ------- ------- -------
 **1X8**    nosep     core     nosep   nosep   nosep   nosep   nosep   nosep 

 **1X10**   nosep    nosep     nosep   nosep   nosep   nosep   nosep   nosep 

 **1X17**   nosep     core     nosep   nosep   nosep   nosep   nosep   nosep 

 **1X37**   nosep   outsider   nosep   nosep   nosep   nosep   nosep   nosep 

 **1X39**   nosep    nosep     nosep   nosep   nosep   nosep   nosep   nosep 

 **1X54**   nosep    nosep     nosep   core    nosep   nosep   nosep   nosep 
-----------------------------------------------------------------------------

Table: Table continues below

 
-------------------------------------------------------------------------
  &nbsp;     X52     X53     X54     X55     X56     X58     X59    X60  
---------- ------- ------- ------- ------- ------- ------- ------- ------
 **1X8**    nosep   nosep   nosep   nosep   nosep   nosep   core    core 

 **1X10**   nosep   nosep   nosep   nosep   nosep   nosep   core    core 

 **1X17**   nosep   nosep   nosep   nosep   nosep   nosep   core    core 

 **1X37**   nosep   nosep   nosep   nosep   nosep   nosep   core    core 

 **1X39**   nosep   nosep   nosep   nosep   nosep   nosep   nosep   core 

 **1X54**   nosep   nosep   nosep   nosep   nosep   nosep   core    core 
-------------------------------------------------------------------------

Table: Table continues below

 
------------------------------------------------------------
  &nbsp;     X62     X65     X66     X67      X68      X69  
---------- ------- ------- ------- ------- ---------- ------
 **1X8**    nosep   core    core    nosep     core     core 

 **1X10**   nosep   core    core    nosep     core     core 

 **1X17**   nosep   core    core    nosep     core     core 

 **1X37**   nosep   nosep   nosep   nosep   outsider   core 

 **1X39**   nosep   core    core    nosep    nosep     core 

 **1X54**   nosep   core    core    nosep     core     core 
------------------------------------------------------------

```r
#create a parallel cluster with the cluster log output
if(logfile %>% file.exists()){
        file.remove(logfile)
}
```

```
## [1] TRUE
```

```r
cluster <<- makeCluster(6,outfile = logfile)
registerDoParallel(cluster)
#create a diriving summary of genes
master.table.core.per.gene.summary<-as.data.frame(
        t(
                parApply(cl=cluster, X = master.table.core, MARGIN = 1, FUN = function(i){
                        sum.core<-length(i[i=="core"])
                        sum.outsider<-length(i[i=="outsider"])
                        sum.nosep<-length(i[i=="nosep"])
                        return(c(core=sum.core, outsider=sum.outsider, nonsep=sum.nosep))
})))
#do not forget to close the cluster
stopCluster(cluster)
```

The so called "frequent outsides" below are essentially the genes that were marked by the models as "outsiders" many times, thus suggesting that the origin of these orfs is different from their host genome.

```r
#append the names because the following filter will erase the rownames
master.table.core.per.gene.summary$gene.name<-row.names(master.table.core.per.gene.summary)
#separate the "always core table"
master.table.core.always.core<- master.table.core.per.gene.summary %>% filter(outsider==0)
head(master.table.core.always.core) %>% pander(caption="Mater Table Core")
```


--------------------------------------
 core   outsider   nonsep   gene.name 
------ ---------- -------- -----------
  17       0         44        1X8    

  12       0         49       1X10    

  17       0         44       1X17    

  8        0         53       1X39    

  13       0         48       1X54    

  22       0         39       1X59    
--------------------------------------

Table: Mater Table Core

```r
#separate the "outsider once"
master.table.core.outsider.once<- master.table.core.per.gene.summary %>% filter(outsider==1)
head(master.table.core.outsider.once) %>% pander(caption="One Time Outsiders")
```


--------------------------------------
 core   outsider   nonsep   gene.name 
------ ---------- -------- -----------
  4        1         56       1X638   

  8        1         52      1X2086   

  3        1         57      1X2917   

  9        1         51      1X7396   

  4        1         56      1X9503   

  7        1         53      1X10411  
--------------------------------------

Table: One Time Outsiders

```r
#separate the frequent outsider
master.table.core.outsider.frequently<- master.table.core.per.gene.summary %>% filter(outsider>1)
head(master.table.core.outsider.frequently) %>% pander(caption="Frequent Outsiders")
```


--------------------------------------
 core   outsider   nonsep   gene.name 
------ ---------- -------- -----------
  6        9         46       1X37    

  3        7         51      1X1773   

  4        6         51      1X2780   

  3        12        46      1X2789   

  6        3         52      1X3324   

  6        7         48      1X3444   
--------------------------------------

Table: Frequent Outsiders

```r
#write out the data
write.table(master.table.core.always.core, 
            file="master.table.core.always.core.bootstrapped.txt" %>% sixtyseven.genomes.quadruple.model, 
            sep="\t", col.names=TRUE, quote = FALSE)
write.table(master.table.core.outsider.once, 
            file="master.table.core.outsider.once.bootstrapped.txt" %>% sixtyseven.genomes.quadruple.model, 
            sep="\t", col.names=TRUE, quote = FALSE)
write.table(master.table.core.outsider.frequently, 
            file="master.table.core.outsider.frequently.bootstrapped.txt" %>% sixtyseven.genomes.quadruple.model, 
            sep="\t", col.names=TRUE, quote = FALSE)
#
```

The hectic superlogng function below is a minor modification of the one that is present in the `gbra` helper package and is redefined here because the previous one was very slow, and this one is faster (and will replace the original one in `gbra`).

```r
calculate.clustering.table<- function(bh.data, rsquared.cutoff, rsquared.table, legend, processors=4, data.folder=data.folder){
        registerDoMC(processors)
        #for each match table find the genes that are in the core and those, that are outsiders
        core.outsider.list<-foreach(x=1:nrow(rsquared.table), .errorhandling = "stop", .verbose = FALSE)%dopar%{
                #extract i and j from the filename
                ij<- strsplit(x = rsquared.table$match.file.name[x], split = ".", fixed=TRUE)
                ij<- strsplit(x = ij[[1]][1], split = "_", fixed=TRUE)
                i<- as.numeric(ij[[1]][1])
                j<- as.numeric(ij[[1]][2])
                if(file.exists(data.folder(rsquared.table$match.file.name[x]))){
                         match.table<- fread(input = data.folder(rsquared.table$match.file.name[x]), sep = "\t", data.table = FALSE, stringsAsFactors = FALSE)
                }else{
                         match.table<- select.genome.data(bh.data, c(i,j)) %>% filter.genome.data(c(i,j))
                         match.table<- match.table %>% mutate(name=rownames(match.table))
                }

                #split by genome ids
                i.columns<- match.table %>% filter(grepl(pattern = paste0("\\b",i,"X"),x = name,perl = TRUE))
                j.columns<- match.table %>% filter(grepl(pattern = paste0("\\b",j,"X"),x = name,perl = TRUE))
                #if the rsquared is too low, then just return a nonsep-filled columns
                if(is.na(rsquared.table$rsquared[x])||rsquared.table$rsquared[x]<rsquared.cutoff){
                        i.col=rep("nosep",nrow(i.columns))
                        names(i.col)<- i.columns$name
                        j.col=rep("nosep",nrow(j.columns))
                        names(j.col)<- j.columns$name
                        return(list(i.col=i.col, j.col=j.col, i=i, j=j))
                }else{
                        i.col=sapply(i.columns$match, function(k){
                                if(is.na(k)){
                                        return("nosep")
                                }
                                if(k){
                                        return("core")
                                }else{
                                        return("outsider")
                                }
                        })
                        names(i.col)<- i.columns$name
                        j.col=sapply(j.columns$match, function(k){
                                if(is.na(k)){
                                        return("nosep")
                                }
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
```

The rest of the funcitons can be borrowed form previous analysis.

```r
#first prepare the functions
source("R/extract.cores.R")
source("R/convert.core.list.matrix.R")
```

The enormous and unreadable pipe below will create a "clustering list" (that is a list, which is essentially an upper triangular matrix of orfs vs genomes).

```r
#load the core-outsider data
model.df<- cbind.data.frame(query.genome=model.driver$X1, target.genome=model.driver$X2, rsquared=model.driver$rsquared)
model.df$match.file.name<- sapply(1:nrow(model.df), function(i){
        return(
                paste0(model.df[i,1], '_', model.df[i,2], '.txt.bootstrapped.match')
        )
})
model.df[is.na(model.df[,3]),3]<-0

clustering.list<- calculate.clustering.table(
  bh.data = bh.data, rsquared.cutoff = 0.73, rsquared.table = model.df, legend = legend, processors = 8,
  data.folder=sixtyseven.genomes.quadruple.model.match)

core.distribution.data<- convert.core.list.matrix(clustering.list = clustering.list, legend = legend)
```

```
## [1] "list"
```

A plot below will demonstate the amount of genes from outsider and core groups. If a model had an R-squared less than 0.73, then all genes are greyed out as "no separation possible".

```r
rsq.bar<- core.distribution.data %>% select.genome.data(6) %>% filter.genome.data(6) %>% mutate(id=rownames(.)) %>% melt(id.vars="id") %>% group_by(variable,value) %>%
  summarise(times=n()) %>% do({
    df<- .
    df$times<- rescale(df$times, from = c(0,sum(df$times)),to = c(0,1))
    df
  }) %>% ungroup() %>% mutate(variable=str_replace_all(string=.$variable, pattern ="X", replacement="")) %>%
  merge(x=., y=legend, by.x="variable", by.y="id_genomes") %>% 
  merge(x=., y=filter(model.df,query.genome==6|target.genome==6) %>% mutate(variable=sapply(1:nrow(.),function(x){
    if(.$target.genome[x]==6){
      return(.$query.genome[x])
    }else{
      return(.$target.genome[x])
    }
  }))) %>% arrange(rsquared) %>% 
  mutate(name=str_trim(toupper(str_replace_all(string=.$name, pattern = "\\.fasta|_|complete_genome|chromosome", replacement=" ")))) %>%
  mutate(name=factor(.$name, levels=unique(.$name))) %>%
  ggplot(data=., mapping=aes(x=name, y=times, group=variable, fill=value, order=value)) +geom_bar(stat="identity") + 
  theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1),plot.margin = unit(c(1, 1, 1, 4), "cm")) +
  xlab("Genome") + ylab("Core/Outsider Ratio") + labs(title="Drosophila Core/Outsider Gene Classification") +
  geom_segment(mapping=aes(x=1, xend=nrow(legend)-1, y=0.1, yend=0.1), arrow=arrow(length=unit(0.5,"cm")),colour="red") +
  scale_fill_manual(values=c("#06960B","#999999","#F7728A"),labels=c("Core Genes","No Clear Separation","Outsider Genes"), guide=guide_legend(title=""))

plot(rsq.bar)
```

![](bootstrapped.models_files/figure-html/core outsider plot for drosophila-1.png)<!-- -->

```r
pdf(file="core.outsider.barplot.bootstrapped.pdf" %>% sixtyseven.genomes.quadruple.model, width=15, height=6)
plot(rsq.bar)
dev.off()
```

```
## png 
##   2
```

Now, the famous plot between *D.simulans* and *Wolbachia* that will show two clusters and brign red "outsiders" (candidates for HGT) from *Wolbachia* (would have been from *D.simulans* if ther were any).

```r
dph.arch.pair<- c(3, 6)
dph.wb<- bh.data %>% select.genome.data(dph.arch.pair) %>% filter.genome.data(dph.arch.pair)
dph.wb.clust.one<- (core.distribution.data %>% 
        select.genome.data(dph.arch.pair, exclude=FALSE) %>% 
        filter.genome.data(dph.arch.pair[1]))[,2]
dph.wb.clust.two<- (core.distribution.data %>% 
        select.genome.data(dph.arch.pair, exclude=FALSE) %>% 
        filter.genome.data(dph.arch.pair[2]))[,1]
dph.wb.clust<- c(dph.wb.clust.one, dph.wb.clust.two)
dph.wb.pca<- prcomp(dph.wb)
dph.wb.pca.plot.data<- cbind.data.frame(GENE=rownames(dph.wb),
                                        ORGANISM=sapply(rownames(dph.wb),function(x){
                                                return(
                                                        strsplit(x = x, split = "X", fixed = TRUE)[[1]][1]
                                                )
                                        }), 
                                        PC1=dph.wb.pca$x[,1], 
                                        PC2=dph.wb.pca$x[,2], 
                                        CLUSTER=dph.wb.clust)
dph.wb.pca.plot.data$STATE<- sapply(1:nrow(dph.wb.pca.plot.data), i:={
    if(dph.wb.pca.plot.data[i,]$CLUSTER == 'nosep'){
        return('Uncertain')
    }
    if(dph.wb.pca.plot.data[i,]$CLUSTER == 'core' && dph.wb.pca.plot.data[i,]$ORGANISM == 3){
      return('Wolbachia Core')
    }
    if(dph.wb.pca.plot.data[i,]$CLUSTER == 'outsider' && dph.wb.pca.plot.data[i,]$ORGANISM == 3){
      return('Wolbachia Outsider')
    }
    if(dph.wb.pca.plot.data[i,]$CLUSTER == 'core' && dph.wb.pca.plot.data[i,]$ORGANISM == 6){
      return('Drosophila Core')
    }
    if(dph.wb.pca.plot.data[i,]$CLUSTER == 'outsider' && dph.wb.pca.plot.data[i,]$ORGANISM == 6){
      return('Drosophila Outsider')
    }
})
dph.wb.pca.plot<- dph.wb.pca.plot.data %>% ggplot(data=., mapping=aes(x=PC1, y=PC2)) +
        geom_point(aes(color=STATE), alpha=0.5) + theme_bw() +
        scale_color_manual(values=c("lightpink", "red", "lightgray", "lightblue", "blue"),
                           guide = guide_legend(title="CLUSTER")) + theme(legend.key = element_blank())
dph.wb.pca.plot.facet<- dph.wb.pca.plot + facet_wrap(~STATE,nrow = 1)
plot(dph.wb.pca.plot)
```

![](bootstrapped.models_files/figure-html/a more detailed drosophila - wolbachia plot-1.png)<!-- -->

```r
pdf(file = "w.eds.vs.ds.bootstrapped.pdf" %>% sixtyseven.genomes.quadruple.model, width=10, height=10)
dev.off()
```

```
## png 
##   2
```

```r
pdf(file = "w.eds.vs.ds.bootstrapped.facet.pdf" %>% sixtyseven.genomes.quadruple.model, width=25, height=5)
plot(dph.wb.pca.plot.facet)
dev.off()
```

```
## png 
##   2
```

And a totally useless table that shows proportions of "core", "outsider", "nosep" for each genome. On the second thought, it is not as useless if one by any chance is interested which genome has very bad separation and thus is not really "very unique" so to speak (meaning looks like all the others in a way by its gene compound).

```r
#extract the genomes
bar.plot.core.outsider<- expand.df(short.df = core.distribution.data) %>% select(-QUERY_ORF_ID) %>% 
  merge(x=., y=legend, by.x="ID_QUERY_GENOME", by.y="id_genomes")%>% select(-ID_QUERY_GENOME) %>%
  melt(id.vars="name") %>% group_by(name,value) %>% summarize(total=n()) %>% ungroup() %>% group_by(name) %>%
  do({
    df<- .
    df$total=rescale(df$total, from=c(0, sum(df$total)), to=c(0,1))
    df
  }) %>% ungroup() %>% mutate(name=sapply(name, function(x){
    if(x=="WOLBACHIA ENDOSYMBIONT OF DROSOPHILA SIMULANS WHA"){
      return("WOLBACHIA DROSOPHILA")
    }
    if(x=="WOLBACHIA ENDOSYMBIONT OF CULEX QUINQUEFASCIATUS PEL"){
      return("WOLBACHIA CULEX")
    }
    return(x)
  })) %>% 
  mutate(name=sapply(name,function(x){
    x<- strsplit(x,split = " ",fixed = TRUE)
    return(paste(x[[1]][1],x[[1]][2],sep=" "))
    })) %>%
  ggplot(data=., mapping=aes(x=value, y=total, group=name, fill=name)) + geom_bar(stat="identity", position="dodge") +
  theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1, size=24), legend.position="none", title= element_text(size=32))+
  facet_wrap(~name,ncol=4) + xlab("") + ylab("proportion") + labs(title="Core/Outsider/NoSep Distribution Over Genomes")
plot(bar.plot.core.outsider)
```

![](bootstrapped.models_files/figure-html/full core distribution plot-1.png)<!-- -->

```r
pdf(file="bar.plot.core.outsider.bootstrapped.pdf" %>% sixtyseven.genomes.quadruple.model, width=20, height=35)
plot(bar.plot.core.outsider)
dev.off()
```

```
## png 
##   2
```
