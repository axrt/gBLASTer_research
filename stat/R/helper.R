#' Use this function to rename the trees from the short codes to an appropriate legend names
rename.tree<-function(tree, legend){
  tree$tip.label<-sapply(tree$tip.label,function(i){
    if(grepl(x=i, pattern = "X",fixed = TRUE)){
      return(paste0(legend$name[which(legend$id_genomes==extract_numeric(i))],"_genome"))
    }else{
      return(paste0(legend$name[which(legend$id_genomes==extract_numeric(i))],"_OUTSIDER>>>>>>"))
    }
  })
  return(tree)
}

tidy.tree<-function(tree, legend){
  tree$tip.label<-sapply(tree$tip.label,function(i){
    if(grepl(x=i, pattern = "X",fixed = TRUE)){
      return(legend$name[which(legend$id_genomes==extract_numeric(i))])
    }else{
      return(legend$name[which(legend$id_genomes==extract_numeric(i))])
    }
  })
  return(tree)
}

#' 
fitch.output.folder<- function(file=""){
  return(paste("fitch.output", file, sep="/"))
}

#'
data.folder<- function(file=""){
  return(paste("data",file,sep="/"))
}

#'
data.bh<- function(file=""){
  return(paste("..","bh",file, sep="/"))
}

#'
sixtyone.genomes.dir<- function(file=""){
  return(paste("61.genomes", file, sep="/"))
}

#'
sixtyone.genomes.data.dir<- function(file=""){
  return(paste(sixtyone.genomes.dir("data"), file, sep="/"))
}

#'
sixtyone.genomes.pred.dir<- function(file=""){
  return(paste(sixtyone.genomes.dir("pred"), file, sep="/"))
}
#'
#'
sixtyseven.genomes.dir<- function(file=""){
  return(paste("67.genomes", file, sep="/"))
}

#'
sixtyseven.genomes.data.dir<- function(file=""){
  return(paste(sixtyseven.genomes.dir("data"), file, sep="/"))
}
#'
sixtyseven.genomes.pred.dir<- function(file=""){
  return(paste(sixtyseven.genomes.dir("pred"), file, sep="/"))
}
#'
sixtyseven.genomes.quadruple.model<- function(file=""){
  return(paste("67.genomes.quadruple.model", file, sep="/"))
}
#'
sixtyseven.genomes.quadruple.model.models<- function(file=""){
  return(paste(sixtyseven.genomes.quadruple.model("models"), file, sep="/"))
}

#'
sixtyseven.genomes.quadruple.model.bootstrapped<- function(file=""){
        return(paste(sixtyseven.genomes.quadruple.model("bootstrapped"), file, sep="/"))
}

#'
sixtyseven.genomes.quadruple.model.match<- function(file=""){
        return(paste(sixtyseven.genomes.quadruple.model("match"), file, sep="/"))
}

#'
tidy.legend<- function(legend){
  legend$name<- str_trim(toupper(str_replace_all(string=legend$name, pattern = "\\.fasta|_|complete_genome|chromosome", replacement=" ")))
  return(legend)
}
#`
research.dir<- function(file=""){
  return(paste("../research",file,sep="/"))
}
#`
research.domains.dir<- function(file=""){
  return(research.dir(paste("domain", file,sep = "/")))
}
#'

