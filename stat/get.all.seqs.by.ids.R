source("R/helper.R")
library("gbra")

load(sixtyseven.genomes.dir("bh.data.normal.rda"))
connection<- connect.derby(db = "/home/alext/Documents/Research/gBLASTer/db/gblasterdb", usr = "APP",
                           derby.jar = "/home/alext/Developer/gBLASTer/out/artifacts/gBLASTer_jar/derby-10.9.1.0.jar")
write.out.seqs<- function(connection, ids, file, sep="\t"){
  for(id in ids){
    sequence<- dbGetQuery(conn = connection, statement = paste0("select sequence from orfs where id_orf=",id))
    write(paste(id,sequence$SEQUENCE,seq=sep), file, append=TRUE)
  }
  return(file)
}

write.out.seqs(connection$conn, expand.df(bh.data.norm)$QUERY_ORF_ID, file=sixtyseven.genomes.dir("bh.data.norm.seqs.txt")) 

load(sixtyseven.genomes.dir("core.distribution.data.rda"))
write.table(core.distribution.data, file = sixtyseven.genomes.dir("core.distribution.data.txt"), sep = "\t", quote = FALSE)
write.table(core.distribution.data.melt.summary.total, file = sixtyseven.genomes.dir("core.distribution.data.summary.total.txt"), sep = "\t", quote = FALSE)
