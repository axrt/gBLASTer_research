# Minimum Hit Cutoff
Alexander Tuzhikov  
September 30, 2017  
Here I would like to take a look at the distribution of the hits.

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
library(reshape2)
```

We create a legend, that will only contain the genomes that we want to include (the others are not good enough for this or that reason)

```r
legend<- read.table(file=sixtyseven.genomes.dir("legend.txt"), sep="\t", header=TRUE, stringsAsFactors = FALSE) %>% 
  filter(id_genomes!=23) %>% 
  filter(id_genomes!=40) %>% 
  filter(id_genomes!=47) %>%
  filter(id_genomes!=57) %>%
  filter(id_genomes!=63) %>%
  filter(id_genomes!=64) %>%
  select(name, id_genomes) %>%
  arrange(id_genomes)
```

We load the bidirectional hits form the previous iterations of analaysis.

```r
load(sixtyseven.genomes.dir(bh.data.raw.file))
```


```r
bh.data %>%
  restrict.minimal.hits(minhit = 2) %>%
  as.data.frame -> bh.data.mask
bh.data.mask[,3:ncol(bh.data.mask)][bh.data.mask[,3:ncol(bh.data.mask)]>0]<-1
bh.data.hits<- data.frame(genome=bh.data.mask$ID_QUERY_GENOME, 
                          number_of_hits = rowSums(bh.data.mask[,3:ncol(bh.data.mask)]),
                          poisson.vals=rpois(n=nrow(bh.data.mask), lambda = 1)+1)
bh.data.hits %>% ggplot(aes(x=number_of_hits)) + geom_histogram(bins=67)
```

![](minhit.cutoff.distribution_files/figure-html/sums of hits-1.png)<!-- -->

```r
bh.data.hits %>% ggplot(aes(x=poisson.vals)) + geom_histogram(bins=67) + xlim(c(1, 10))
```

```
## Warning: Removed 1 rows containing missing values (geom_bar).
```

![](minhit.cutoff.distribution_files/figure-html/sums of hits-2.png)<!-- -->

```r
bh.data.hits %>% group_by(poisson.vals) %>% summarize(poisson=n()) -> poisson.summary
bh.data.hits %>% group_by(number_of_hits) %>% summarize(gblaster=n()) -> numhits.summary
cbind.data.frame(rbind.data.frame(poisson.summary[-1,], data.frame(poisson.vals=10:61, poisson = 0)), numhits.summary)[,-1] %>%
  reshape2::melt(id.vars='number_of_hits') %>% ggplot(aes(x=number_of_hits, y=value, group=variable, color=variable)) + geom_line() + scale_x_continuous(breaks=seq(0, 67,2))
```

![](minhit.cutoff.distribution_files/figure-html/sums of hits-3.png)<!-- -->

```r
bh.data.hits %>% group_by(genome) %>% do({
  df<-.
  df$number_of_hits<- df$number_of_hits/nrow(df)
  df
}) -> bh.data.hits.normalized.by.genome.size
bh.data.hits.normalized.by.genome.size %>% ggplot(aes(x=number_of_hits)) + geom_histogram(bins=67)
```

![](minhit.cutoff.distribution_files/figure-html/sums of hits-4.png)<!-- -->

```r
mean(bh.data.hits$number_of_hits)
```

```
## [1] 7.686517
```

```r
sd(bh.data.hits$number_of_hits)
```

```
## [1] 10.15242
```
