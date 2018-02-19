total.logit.model<- function(df,
                             legend,cores = 1,
                             sep = "X", 
                             outfile = "cluster.log", 
                             do.step=TRUE, 
                             direction="backward") {
    if (!require("doParallel")) {
      install.packages("doParallel")
      library("doParallel")
    }
    if (!require("mgcv")) {
      install.packages("mgcv")
      library("mgcv")
    }
    
    env <- environment()
    
    #calculate the number of models that will be trained
    number.to.run <- upper.triangle.diag(length = length(legend$id_genomes))
    message("Logregs to run:", number.to.run)
    
    #initialize a cluater that will carry out the training
    cluster <<- makeCluster(cores, outfile = outfile)
    registerDoParallel(cluster)
    
    out.list <- list()
    comb <- function(x,y) {
      x[[length(x)+1]]<-y
      return(x)
    }
    dir<-"data"
    if(!file.exists(dir)){
      dir.create("data")
    }
    for (i in 1:(nrow(legend) - 1)) {
      out.list[[i]]<-foreach(
        j = (i + 1):nrow(legend),.combine = "list",
        .inorder = TRUE,.packages = "gbra",.errorhandling = "stop"
      ) %dopar% {
        out.file<-paste(dir,paste(paste(legend$id_genomes[i],legend$id_genomes[j],sep="_"),".rda",sep = ""),sep = "/")
        if(!file.exists(out.file)){
            l.m.g<-logreg.mismatch.genes(
            data = df, org1 = legend$id_genomes[i], org2 = legend$id_genomes[j], keep.mod=TRUE, do.step=do.step, direction=direction
            )
            save(l.m.g,file = out.file)
            return(list(glm.summary=summary(l.m.g),Rsq=l.m.g[5]))
        }else{
          return(NULL)
        }
      }
    }
    
    stopCluster(cluster)
    return(out.list)
  }

#'This is a rather dumb function to calculate the number of members of the upper triangle matrix + its diagonal.
#'Used to calculate the number of genome pairs to have models trained on.
#'@param length positive integer
#'@return integer, number of members in the upper triangle + diagonal
upper.triangle.diag<- function(length){
  return(
    (length^2) - (length/2) 
  )
}