for(block in 2:40){
  betaML_file <- paste0("simulation_results/coverage/beta_ML_coverage_block",block,".csv")
  beta <- read.csv(betaML_file,header=F)
  beta[(26*432+1):nrow(beta),1] <- beta[(26*432+1):nrow(beta),1] + 1
  write.table(beta,file=betaML_file,row.names=F,col.names=F,sep=",")
  
  VI_file <- paste0("simulation_results/coverage/VI_coverage_block",block,".csv")
  VI <- read.csv(VI_file,header=F)
  VI[(26*571+1):nrow(VI),1] <- VI[(26*571+1):nrow(VI),1] + 1
  write.table(VI,file=VI_file,row.names=F,col.names=F,sep=",")
  
  hyper_file <- paste0("simulation_results/coverage/hyper_coverage_block",block,".csv")
  hyper <- read.csv(hyper_file,header=F)
  hyper[27:nrow(hyper),1] <- hyper[27:nrow(hyper),1] + 1
  write.table(hyper,file=hyper_file,row.names=F,col.names=F,sep=",")
  
}