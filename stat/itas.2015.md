# ITAS 2015
Alexander Tuzhikov  
September 4, 2015  



```r
#yr libraries
source("R/yr.lib.R")
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
## The following objects are masked from 'package:stats':
## 
##     filter, lag
## 
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
## 
## Loading required package: tidyr
```

```r
#hmptrees
source("R/hmptrees.lib.R")
```

```
## Loading required package: HMPTrees
## Loading required package: ape
```

```r
#and the extension
source("R/hmptreesextend.lib.R")
```

```
## Loading required package: devtools
## Loading required package: HMPTreesExtend
```

```r
#rphilyp
source("R/rphylip.lib.R")
```

```
## Loading required package: Rphylip
```

```r
#gbra
library("gbra")
```

```
## Loading required package: data.table
## 
## Attaching package: 'data.table'
## 
## The following objects are masked from 'package:dplyr':
## 
##     between, last
## 
## Loading required package: Rcpp
## Loading required package: boot
## Loading required package: stringr
```

```r
#ggplot
source("R/ggplot2.lib.R")
```

```
## Loading required package: ggplot2
```

```r
#reshape2
source("R/reshape2.lib.R")
```

```
## Loading required package: reshape2
```

```r
#httr
source("R/httr.lib.R")
```

```
## Loading required package: httr
```

```r
#stringr
source("R/stringr.lib.R")
#grid
source("R/gridextra.lib.R")
```

```
## Loading required package: gridExtra
## Loading required package: grid
```

```r
#scales
source("R/scales.lib.R")
```

```
## Loading required package: scales
```

```r
#doMC
source("R/domc.lib.R")
```

```
## Loading required package: doMC
## Loading required package: foreach
## Loading required package: iterators
## Loading required package: parallel
```

```r
#doparallel
source("R/doparallel.lib.R")
```

```
## Loading required package: doParallel
```

```r
#data.table
source("R/data.table.lib.R")
```


```r
source("R/helper.R")
#bh data with 4 cutoff
bh.data<- read.table(file="bh.data.minhit_4.txt")#data
legend<-tidy.legend(load.legend(legend.file = "/home/alext/Developer/gbra/inst/extradata/legend.csv",sep = "\t")) #legend
model.df<- read.table(file = "model.df.txt", header = TRUE, sep = "\t")#model descriptions
#folder helper
self.dir<- "itas.2015.output"
if(!file.exists(self.dir)){
  dir.create(self.dir)
}
self.dir.folder<- function(file=""){
  return(paste("itas.2015.output",file, sep="/"))
}
```


```r
#extract mles from the genomes
dirty.mles<- extract.gprimes(mle.list = getMLEs(df = expand.df(bh.data)))
```

```
## Calculating MLE for genome 1
## Calculating MLE for genome 2
## Calculating MLE for genome 3
## Calculating MLE for genome 4
## Calculating MLE for genome 5
## Calculating MLE for genome 6
## Calculating MLE for genome 7
## Calculating MLE for genome 8
## Calculating MLE for genome 9
## Calculating MLE for genome 10
## Calculating MLE for genome 11
## Calculating MLE for genome 12
## Calculating MLE for genome 13
## Calculating MLE for genome 14
## Calculating MLE for genome 15
## Calculating MLE for genome 16
## Calculating MLE for genome 17
## Calculating MLE for genome 18
## Calculating MLE for genome 19
## Calculating MLE for genome 20
## Calculating MLE for genome 21
## Calculating MLE for genome 22
## Calculating MLE for genome 23
## Calculating MLE for genome 24
## Calculating MLE for genome 25
## Calculating MLE for genome 26
## Calculating MLE for genome 27
## Calculating MLE for genome 28
## Calculating MLE for genome 30
## Calculating MLE for genome 31
## Calculating MLE for genome 32
## Calculating MLE for genome 33
## Calculating MLE for genome 34
## Calculating MLE for genome 35
## Calculating MLE for genome 36
## Calculating MLE for genome 37
## Calculating MLE for genome 38
## Calculating MLE for genome 39
## Calculating MLE for genome 40
## Calculating MLE for genome 41
## Calculating MLE for genome 42
## Calculating MLE for genome 43
## Calculating MLE for genome 44
## Calculating MLE for genome 45
## Calculating MLE for genome 46
```

```r
#run fitch
rf<- Rfitch(D= dist(dirty.mles), path = "/usr/local/bin", model="FM", power=2, negative=TRUE, global=FALSE,random.order=FALSE, root=FALSE, quiet=TRUE)
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
##   !    !  !                +---4         
##   !    !  +----------------2 
##   !  +-3                   +----3         
##   !  ! ! 
##   !  ! !      +-------------34        
##   !  ! +-----32  
##   !  !        +-------------33        
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
## Sum of squares =     0.47957
## 
## Average percent standard deviation =     1.55709
## 
## Between        And            Length
## -------        ---            ------
##    1          2                 0.56870
##    1            31              0.02331
##   31             3              0.02188
##    3            26              0.10007
##   26          28                0.53169
##   26             2              0.57626
##    2          4                 0.13306
##    2          3                 0.14135
##    3            32              0.23553
##   32          34                0.45110
##   32          33                0.44843
##   31             4              0.00843
##    4            30              0.01945
##   30          32                0.61142
##   30            33              0.00654
##   33          35                0.61757
##   33            14              0.04554
##   14            42              0.06021
##   42          44                0.60431
##   42            17              0.01705
##   17            22              0.08082
##   22          24                0.53646
##   22          19                0.54216
##   17            16              0.01402
##   16          18                0.58605
##   16            23              0.07482
##   23          25                0.53224
##   23          16                0.51465
##   14             8              0.06890
##    8            35              0.91406
##   35          37                0.00000
##   35          23                0.00000
##    8             5              0.01146
##    5             9              0.01280
##    9            11              0.00862
##   11            38              0.37322
##   38          40                0.40986
##   38            43              0.03417
##   43          45                0.36484
##   43          7                 0.38096
##   11            25              0.01731
##   25          27                0.63419
##   25            24              0.00742
##   24          26                0.64876
##   24            37              0.00530
##   37          39                0.72057
##   37            20              0.05112
##   20          22                0.63436
##   20             7              0.01657
##    7          9                 0.60368
##    7            15              0.02287
##   15            27              0.01972
##   27          29                0.62623
##   27          17                0.59833
##   15            12              0.06901
##   12            19              0.02489
##   19          21                0.52195
##   19            41              0.06310
##   41          43                0.42543
##   41          14                0.44555
##   12          6                 0.60750
##    9            21              0.00296
##   21            36              0.03480
##   36          38                0.63502
##   36            40              0.01192
##   40          42                0.61509
##   40          11                0.61965
##   21            39              0.00610
##   39          41                0.63752
##   39            13              0.00892
##   13            34              0.11244
##   34          36                0.56428
##   34          15                0.56234
##   13            29              0.03081
##   29          31                0.60576
##   29            18              0.03494
##   18            28              0.00795
##   28          30                0.57362
##   28          20                0.61729
##   18            10              0.03079
##   10          12                0.55188
##   10          10                0.54946
##    5             6              0.01886
##    6          13                0.62382
##    6          8                 0.61395
##    4          5                 0.59296
##    1          1                 0.57093
```

```r
#rename the tree
rf<- tidy.tree(tree = rf, legend = legend)
#save
drity.tree.newick<- "dirty.tree.newick"
write.tree(phy = rf, file = self.dir.folder(drity.tree.newick))
#plot and save results
plot.itol(newick.file = self.dir.folder(drity.tree.newick), output.file = self.dir.folder(paste0(drity.tree.newick,".pdf")))
```

```
## [1] "1291711501384786614415972350"
```

```
## [1] "itas.2015.output/dirty.tree.newick.pdf"
```

```r
#another one with eps, cuz i need one editable in illustrator (pdf is not in the case of iTOL)
plot.itol(newick.file = self.dir.folder(drity.tree.newick), output.file = self.dir.folder(paste0(drity.tree.newick,".eps")),format = "eps", font.size = 120)
```

```
## [1] "1291711501384788014415972370"
```

```
## [1] "itas.2015.output/dirty.tree.newick.eps"
```


```r
#split names into query and target genomes
model.df.genomes<- strsplit(x = str_replace_all(string = model.df$file.name, pattern = "\\.rda", replacement = ""), split = "_", fixed = TRUE)
model.df$query.genome<- as.numeric(unlist(sapply(model.df.genomes,function(x){return(x[[1]])})))
model.df$target.genome<- as.numeric(unlist(sapply(model.df.genomes,function(x){return(x[[2]])})))
model.df$match.file.name<- str_replace_all(string = model.df$file.name, pattern = "rda", replacement = "txt.match")
```

Drosophila is **number 6** in the legend, so let's see for the rsquared distribution for this particular genome.


```r
rsq.plot<- model.df %>% filter(query.genome==6|target.genome==6) %>% mutate(target.genome=sapply(1:length(.$target.genome), function(x){
  if(.$target.genome[x]==6){
    return(.$query.genome[x])
  }else{
    return(.$target.genome[x])
  }
})) %>% select(rsquared, target.genome) %>% merge(x=., y=legend, by.x="target.genome", by.y="id_genomes") %>%
  mutate(name=str_trim(toupper(str_replace_all(string=.$name, pattern = "\\.fasta|_|complete_genome|chromosome", replacement=" ")))) %>% 
  arrange(rsquared) %>% mutate(name=factor(.$name, levels=.$name)) %>% melt(id.vars="name") %>%
  filter(variable%in%c("rsquared")) %>% ggplot(data=., mapping=aes(x=name, y=value, group=variable, color=variable))+geom_point()+geom_smooth()+
  theme_bw()+theme(axis.text.x=element_text(angle = 45, hjust = 1),legend.position="none", plot.margin = unit(c(1, 1, 3, 2), "cm"), axis.title.x=element_text(vjust=1)) + 
  geom_hline(yintercept = 0.73, color="red") + 
  xlab("Genome") + ylab("R-squared") + labs(title = "Drosophila simulans models R-squared values")
plot(rsq.plot)
```

```
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
```

![](itas.2015_files/figure-html/rsquared distribution-1.png) 

```r
pdf(file = self.dir.folder("rsq.drosoph.pdf"), width = 15, height=8)
plot(rsq.plot)
```

```
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
```

```r
dev.off()
```

```
## png 
##   2
```


```r
#first prepare the functions
source("R/calculate.clustering.table.R")
source("R/convert.core.list.matrix.R")
source("R/extract.cores.R")
#load the core-outsider data
core.distribution.data<- convert.core.list.matrix(clustering.list = calculate.clustering.table(
  bh.data = bh.data, rsquared.cutoff = 0.73, rsquared.table = model.df, legend = legend, processors = 4), legend = legend)
```

```
## [1] "list"
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
## Warning in rbind_all(clustering.list.cbind): Unequal factor levels:
## coercing to character
```

```r
rsq.bar<- core.distribution.data %>% select.genomes(df=., g.ids=6) %>% mutate(id=rownames(.)) %>% melt(id.vars="id") %>% group_by(variable,value) %>%
  summarise(times=n()) %>% do({
    .$times<- rescale(.$times, from = c(0,sum(.$times)),to = c(0,1))
    return(.)
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

![](itas.2015_files/figure-html/drosophila clusters-1.png) 

```r
pdf(file=self.dir.folder("core.outsider.barplot.pdf"),width=15, height=6)
plot(rsq.bar)
dev.off()
```

```
## png 
##   2
```

Now we are going to use part of the above code to plot out an idealized tree.


```r
#select only core genes
bh.data.cores<- bh.data[extract.cores(bh.data = bh.data, core.number.cutoff=5, core.outsider.df = core.distribution.data, processors = 6),]
#extract mle gprimes
bh.data.cores.mles.gprimes<- extract.gprimes(mle.list = getMLEs(df = expand.df(bh.data.cores)))
```

```
## Calculating MLE for genome 1
## Calculating MLE for genome 2
## Calculating MLE for genome 3
## Calculating MLE for genome 4
## Calculating MLE for genome 5
## Calculating MLE for genome 6
## Calculating MLE for genome 7
## Calculating MLE for genome 8
## Calculating MLE for genome 9
## Calculating MLE for genome 10
## Calculating MLE for genome 11
## Calculating MLE for genome 12
## Calculating MLE for genome 13
## Calculating MLE for genome 14
## Calculating MLE for genome 15
## Calculating MLE for genome 16
## Calculating MLE for genome 17
## Calculating MLE for genome 18
## Calculating MLE for genome 19
## Calculating MLE for genome 20
## Calculating MLE for genome 21
## Calculating MLE for genome 22
## Calculating MLE for genome 23
## Calculating MLE for genome 24
## Calculating MLE for genome 25
## Calculating MLE for genome 26
## Calculating MLE for genome 27
## Calculating MLE for genome 28
## Calculating MLE for genome 30
## Calculating MLE for genome 31
## Calculating MLE for genome 32
## Calculating MLE for genome 33
## Calculating MLE for genome 34
## Calculating MLE for genome 35
## Calculating MLE for genome 36
## Calculating MLE for genome 37
## Calculating MLE for genome 38
## Calculating MLE for genome 39
## Calculating MLE for genome 40
## Calculating MLE for genome 41
## Calculating MLE for genome 42
## Calculating MLE for genome 43
## Calculating MLE for genome 44
## Calculating MLE for genome 45
## Calculating MLE for genome 46
```

```r
#distance
bh.data.cores.mles.gprimes.dist<- dist(bh.data.cores.mles.gprimes)
#fitch
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
```

```r
#normalize the genome names
rf.names<- tidy.tree(tree = rf, legend = legend)
output.file<- self.dir.folder("clean.tree.0.73.rsq.5.corecut.newick")
write.tree(phy = rf.names, file = output.file)
plot.itol(newick.file = output.file, output.file = paste0(output.file,".pdf"))
```

```
## [1] "1291711501384817214415973580"
```

```
## [1] "itas.2015.output/clean.tree.0.73.rsq.5.corecut.newick.pdf"
```

```r
plot.itol(newick.file = output.file, output.file = paste0(output.file,".eps"), format = "eps", font.size = 120)
```

```
## [1] "1291711501384818614415973600"
```

```
## [1] "itas.2015.output/clean.tree.0.73.rsq.5.corecut.newick.eps"
```




















