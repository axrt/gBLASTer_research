master.table.8.3x6.gen<-attach.genome(master.table.8.3x6)
library(ggplot2)
library(reshape2)

master.table.8.3x6.gen.melt<-melt(data = master.table.8.3x6.gen,id.vars = c("genome","names"))
multi.gen.plot<-ggplot(data = master.table.8.3x6.gen.melt, mapping = aes(x=variable, y=value, colour=variable, group=genome))
pdf(file = "3.6.mean.score.pdf",width = 15,height = 15)
multi.gen.plot+stat_summary(aes(group=genome, fill=variable), fun.y=mean, geom="bar")
dev.off()
