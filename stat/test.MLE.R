test.data.frame<-master.table.8[1:100,]
result.data<-data.frame(mean=rowMeans(test.data.frame),mle=getMLE(df = test.data.frame))
