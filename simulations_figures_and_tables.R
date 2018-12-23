# --------------------------------------------------------------------------------------------- #
# ----------------------- Producing figures and tables for simulations ------------------------ #
# --------------------------------------------------------------------------------------------- #

rm(list=ls())
setwd("/Users/emt380/Documents/PhD_Papers/Air_pollution/R_code/MORETreeS/")

### load test tree ###
load("./moretrees/simulation_inputs/inputs.Rdata")
params$nsamp <- as.character(params$nsamp)
nsims <- 1000

require(clusteval)
require(mclust)
require(igraph)

# Functions needed
source("./moretrees/processing_functions.R")

############ Stitching results together ###############

path_file <- "./Results/Simulations/Dec21_sims/"
params$nsamp2 <- as.character(params$nsamp)
params$nsamp2[params$nsamp==1E05] <- "1e\\+05"
params$nsamp2[params$nsamp==1E06] <- "1e\\+06"

for(i in c(4,8)){ # only need to do this for 4th and 8th sim
  # params
  nsamp <- params$nsamp[i]
  nsamp2 <- params$nsamp2[i]
  whichbeta <- params$whichbeta[i]
  evenness <- params$evenness[i]
  sim_params <- paste("n",nsamp,"_beta",whichbeta,"_evenness",evenness,sep="")
  sim_params2 <- paste("n",nsamp2,"_beta",whichbeta,"_evenness",evenness,sep="")

  # betas
  beta_files <- list.files(path=path_file,pattern=paste("part_betasims_",sim_params2,"_*",sep=""))
  outfile_beta <- paste0(path_file,"betasims_all_",sim_params,".csv")
  file.create(outfile_beta)
  write.table(rbind(c("sim","moretrees_est","uncollapse","truth","adhoc1","adhoc2","adhoc3","adhoc4")), file = outfile_beta, row.names =FALSE, col.names = FALSE,sep = ",", append = TRUE)
  sims_idx <- character(0)
  for(fn in beta_files){
    sims <- read.csv(file=paste0(path_file,fn),header=F,row.names=NULL)
    process_id <- sub(pattern=paste0("part_betasims_",sim_params2,"_"),replace="",x=fn)
    process_id <- paste0("_",sub(pattern=".csv",replace="",x=process_id))
    sims[,1] <- sapply(sims[,1],paste0,process_id)
    sims_idx <- c(sims_idx,unique(sims[,1]))
    write.table(sims, file=outfile_beta, row.names =FALSE, col.names = FALSE,sep = ",", append = TRUE)
  }
  # Keep only first 1000
  sims <- read.csv(file=outfile_beta)
  if(length(sims_idx) > nsims){
    sims <- subset(sims,sim %in% sims_idx[1:nsims])
  }
  sim_n <- unique(sims$sim)
  sim_n <- data.frame(sim=as.character(sim_n),sim_n=1:length(sim_n))
  sims <- merge(sims,sim_n,by="sim",all.x=T)
  outfile_beta <- paste0(path_file,"betasims_",sim_params,".csv")
  file.create(outfile_beta)
  write.table(rbind(c("sim","moretrees_est","uncollapse","truth","adhoc1","adhoc2","adhoc3","adhoc4","sim_n")), file = outfile_beta, row.names =FALSE, col.names = FALSE,sep = ",", append = TRUE)
  write.table(sims, file=outfile_beta, row.names =FALSE, col.names = FALSE,sep = ",", append = TRUE)

  # VI params
  VI_files <- list.files(path=path_file,pattern=paste("part_VIsims_",sim_params2,"_*",sep=""))
  outfile_VI <- paste0(path_file,"VIsims_all_",sim_params,".csv")
  file.create(outfile_VI)
  write.table(rbind(c("sim","mu_gamma","sigma2_gamma","u_s")), file = outfile_VI, row.names=FALSE, col.names=FALSE, sep=",",append=TRUE)
  sims_idx <- numeric(0)
  for(fn in VI_files){
    sims <- read.csv(file=paste0(path_file,fn),header=F,row.names=NULL)
    process_id <- sub(pattern=paste0("part_VIsims_",sim_params2,"_"),replace="",x=fn)
    process_id <- paste0("_",sub(pattern=".csv",replace="",x=process_id))
    sims[,1] <- sapply(sims[,1],paste0,process_id)
    sims_idx <- c(sims_idx,unique(sims[,1]))
    write.table(sims, file=outfile_VI, row.names =FALSE, col.names = FALSE,sep = ",", append = TRUE)
  }
  # Keep only first 1000
  sims <- read.csv(file=outfile_VI)
  sims <- subset(sims,sim %in% sim_n$sim)
  sims <- merge(sims,sim_n,by="sim",all.x=T)
  #if(length(unique(sims$sim)) != nsims) stop(paste("Error in Simulation",i,sep=" "))
  outfile_VI <- paste0(path_file,"VIsims_",sim_params,".csv")
  file.create(outfile_VI)
  write.table(rbind(c("sim","mu_gamma","sigma2_gamma","u_s","sim_n")), file = outfile_VI, row.names=FALSE, col.names=FALSE, sep=",",append=TRUE)
  write.table(sims, file=outfile_VI, row.names =FALSE, col.names = FALSE,sep = ",", append = TRUE)

  # hyperparams
  hyper_files <- list.files(path=path_file,pattern=paste("part_hyperparams_",sim_params2,"_*",sep=""))
  outfile_hyper <- paste0(path_file,"hyperparams_all_",sim_params,".csv")
  file.create(outfile_hyper)
  write.table(rbind(c("sim","rho","tau")), file = outfile_hyper, row.names=FALSE, col.names=FALSE, sep=",",append=TRUE)
  sims_idx <- numeric(0)
  for(fn in hyper_files){
    sims <- read.csv(file=paste0(path_file,fn),header=F,row.names=NULL)
    process_id <- sub(pattern=paste0("part_hyperparams_",sim_params2,"_"),replace="",x=fn)
    process_id <- paste0("_",sub(pattern=".csv",replace="",x=process_id))
    sims[,1] <- sapply(sims[,1],paste0,process_id)
    sims_idx <- c(sims_idx,unique(sims[,1]))
    write.table(sims, file=outfile_hyper, row.names =FALSE, col.names = FALSE,sep = ",", append = TRUE)
  }
  # Keep only first 1000
  sims <- read.csv(file=outfile_hyper)
  sims <- subset(sims,sim %in% sim_n$sim)
  sims <- merge(sims,sim_n,by="sim",all.x=T)
  #if(length(unique(sims$sim)) != nsims) stop(paste("Error in Simulation",i,sep=" "))
  outfile_hyper <- paste0(path_file,"hyperparams_",sim_params,".csv")
  file.create(outfile_hyper)
  write.table(rbind(c("sim","rho","tau","sim_n")), file = outfile_hyper, row.names=FALSE, col.names=FALSE, sep=",",append=TRUE)
  write.table(sims, file=outfile_hyper, row.names =FALSE, col.names = FALSE,sep = ",", append = TRUE)

  # reached_max
  reached_max_files <- list.files(path=path_file,pattern=paste("part_reached_max_",sim_params2,"_*",sep=""))
  outfile_reached_max <- paste0(path_file,"reached_max_all_",sim_params,".csv")
  file.create(outfile_reached_max)
  write.table(rbind(c("sim","reached_max")), file = outfile_reached_max, row.names=FALSE, col.names=FALSE, sep=",",append=TRUE)
  sims_idx <- numeric(0)
  for(fn in reached_max_files){
    sims <- read.csv(file=paste0(path_file,fn),header=F,row.names=NULL)
    process_id <- sub(pattern=paste0("part_reached_max_",sim_params2,"_"),replace="",x=fn)
    process_id <- paste0("_",sub(pattern=".csv",replace="",x=process_id))
    sims[,1] <- sapply(sims[,1],paste0,process_id)
    sims_idx <- c(sims_idx,unique(sims[,1]))
    write.table(sims, file=outfile_reached_max, row.names =FALSE, col.names = FALSE,sep = ",", append = TRUE)
  }
  # Keep only first 1000
  sims <- read.csv(file=outfile_reached_max)
  sims <- subset(sims,sim %in% sim_n$sim)
  sims <- merge(sims,sim_n,by="sim",all.x=T)
  #if(length(unique(sims$sim)) != nsims) stop(paste("Error in Simulation",i,sep=" "))
  outfile_reached_max <- paste0(path_file,"reached_max_",sim_params,".csv")
  file.create(outfile_reached_max)
  write.table(rbind(c("sim","reached_max","sim_n")), file = outfile_reached_max, row.names=FALSE, col.names=FALSE, sep=",",append=TRUE)
  write.table(sims, file=outfile_reached_max, row.names =FALSE, col.names = FALSE,sep = ",", append = TRUE)
  
  # keep track
  print(paste0("Simulation ",i," done."))
}


############ Processing results ########

load("./moretrees/simulation_inputs/inputs.Rdata")
params$nsamp2 <- as.character(params$nsamp)
params$nsamp2[params$nsamp==1E05] <- "1e\\+05"
params$nsamp2[params$nsamp==1E06] <- "1e\\+06"
simres <- data.frame(matrix(nrow=nrow(params),ncol=96))

for(i in (1:nrow(params))){
  # sim parameters
  nsamp <- params$nsamp[i]
  nsamp2 <- params$nsamp2[i]
  whichbeta <- params$whichbeta[i]
  evenness <- params$evenness[i]
  if(whichbeta==1){
    beta_true <- V(tree)$beta1[V(tree)$leaf]
  } else {
    beta_true <- V(tree)$beta2[V(tree)$leaf]
  }
  true.clust <- as.numeric(as.factor(beta_true))
  # load some things
  load(file=paste0("./Data/simdat_n",nsamp2,"_evenness",evenness,".Rdata"))
  # check whether sims converged
  outfile_reached_max <- paste0(path_file,"reached_max_n",nsamp,"_beta",whichbeta,"_evenness",evenness,".csv")
  reached_max <- read.csv(file=outfile_reached_max,header=T,row.names=NULL)
  # extract results
  outfile_beta <- paste0(path_file,"betasims_n",nsamp,"_beta",whichbeta,"_evenness",evenness,".csv")
  betasims <- read.csv(file=outfile_beta,header=T,row.names=NULL)
  nsims.i <- length(unique(betasims$sim))
  outfile_VI <- paste0(path_file,"VIsims_n",nsamp,"_beta",whichbeta,"_evenness",evenness,".csv")
  VIsims <- read.csv(file=outfile_VI,header=T,row.names=NULL)
  print(paste("processing simulation",i,"; number of NAs =", sum(is.na(betasims)),"; number of sims = ",length(unique(betasims$sim)),"; number not converged = ",sum(reached_max$reached_max)))
  # Compute RMSE of grouped estimate for individual betas (ssMOReTreeS)
  group.indiv.rmse <- apply(betasims[,2:8],2,function(x) tapply(x,betasims$sim,rmse.fun,beta=beta_true))
  group.indiv.rmse_tiles <- apply(group.indiv.rmse,2,quantile,c(0.25,0.5,0.75),na.rm=T)
  # Compute RMSE and coverage of individual estimate for individual betas (ssMOReTreeS)
  indiv_betas <- tapply(1:nrow(VIsims),VIsims$sim,indiv.beta.sd.calc,VIsims=VIsims,ancestors=ancestors,pL=pL,p=p,simplify=TRUE)
  indiv_betas_df <- data.frame(beta_indiv=numeric(length=pL*nsims.i),sd_indiv=numeric(length=pL*nsims.i),sim=numeric(length=pL*nsims.i))
  for(j in 1:length(indiv_betas)){
    indiv_betas_df[(pL*(j-1)+1):(pL*j),] <- indiv_betas[j][[1]]
  }
  indiv.rmse.vec <- tapply(indiv_betas_df$beta_indiv,indiv_betas_df$sim,rmse.fun,beta=beta_true)
  indiv.rmse_tiles <- quantile(indiv.rmse.vec,c(0.25,0.5,0.75))
  # Compute number of clusters
  nclust <- tapply(betasims$moretrees_est,betasims$sim,nunique.fun)
  nclust.true <- nunique.fun(true.clust)
  nclust_tiles <- quantile(nclust,c(0.25,0.5,0.75),na.rm=T)
  # Compute ARI 
  ARI <- tapply(betasims$moretrees_est,betasims$sim,function(x) adjustedRandIndex(as.numeric(as.factor(x)),true.clust))
  ARI_tiles <- quantile(ARI,c(0.25,0.5,0.75),na.rm=T)
  ARI_adhoc <- numeric(length=length(groups)*3)
  nclust_adhoc <- numeric(length=length(groups)*3)
  for(j in 1:ncol(groups)){
    g <- groups[,j]
    nclust_adhoc[((j-1)*3+1):(j*3)] <- nunique.fun(g)
    ARI_adhoc[((j-1)*3+1):(j*3)] <- adjustedRandIndex(g,true.clust)
    # added the +1 for the uncollapsed group
  }
  nclust_adhoc <- c(rep(pL,3),nclust_adhoc,rep(pL,3))
  ARI_adhoc <- c(rep(0,3),ARI_adhoc,rep(0,3))
  # Bias of 20 largest coefficients
  bias50 <- apply(cbind(betasims[,2:8],indiv_betas_df$beta_indiv),2,function(x) tapply(X=x,INDEX=betasims$sim,FUN=bias_n_fun,beta_true=beta_true,n=50))
  bias50_tiles <- apply(bias50,2,quantile,p=c(0.25,0.5,0.75),na.rm=T)
  # Put it all into the results data frame
  simres[i,] <- c(group.indiv.rmse_tiles,indiv.rmse_tiles,ARI_tiles,ARI_adhoc,nclust_tiles,nclust_adhoc,bias50_tiles)
}

n.ests <- 8
rmse_names <- sapply(1:(3*n.ests),function(i,x) paste0(x[i,1],x[i,2]),cbind(rep(c("fq_rmse","sq_rmse","tq_rmse"),times=n.ests),as.character(rep(1:n.ests,each=3))))
ARI_names <- sapply(1:(3*n.ests),function(i,x) paste0(x[i,1],x[i,2]),cbind(rep(c("fq_ARI","sq_ARI","tq_ARI"),times=n.ests),as.character(rep(1:n.ests,each=3))))
nclust_names <- sapply(1:(3*n.ests),function(i,x) paste0(x[i,1],x[i,2]),cbind(rep(c("fq_nclust","sq_nclust","tq_nclust"),times=n.ests),as.character(rep(1:n.ests,each=3))))
bias_names <- sapply(1:(3*n.ests),function(i,x) paste0(x[i,1],x[i,2]),cbind(rep(c("fq_bias","sq_bias","tq_bias"),times=n.ests),as.character(rep(1:n.ests,each=3))))
names(simres) <- c(rmse_names,ARI_names,nclust_names,bias_names)

# fq = first quartile of rmse
# sq = second quartile of rmse
# tq = third quartile of rmse

# est1 = moretrees spike & slab
# est2 = uncollapse
# est3 = truth
# est4 = adhoc1
# est5 = adhoc2
# est6 = adhoc3
# est7 = adhoc4 (fully collpsed)
# est8 = moretrees individual estimates

save(simres,params,file="./moretrees/simulation_results/simulation_results_df.Rdata")

######### Creating tables and figures ########

rm(list=ls())
setwd("/Users/emt380/Documents/PhD_Papers/Air_pollution/R_code/MORETreeS/")
load("./moretrees/simulation_inputs/inputs.Rdata")
load("./moretrees/simulation_results/simulation_results_df.Rdata")
simres <- cbind(params,simres)
simres$whichbeta[simres$whichbeta==1] <- "Scenario 1"
simres$whichbeta[simres$whichbeta==2] <- "Scenario 2"
# simres$evenness[simres$evenness==1] <- "even"
# simres$evenness[simres$evenness==2] <- "uneven"
simres$nsamp_lab <- "n=1000"
simres$nsamp_lab[simres$nsamp==10000] <- "n=10,000"
simres$nsamp_lab[simres$nsamp==100000] <- "n=100,000"
simres$nsamp_lab[simres$nsamp==1000000] <- "n=1,000,000"
simres$nsamp_lab <- factor(simres$nsamp_lab,levels=c("n=1000","n=10,000","n=100,000","n=1,000,000"))

require(reshape2)
require(ggplot2)

############ Table A1 - bias of largest estimates ###########

n.ests <- 8
bias_names <- sapply(1:(3*n.ests),function(i,x) paste0(x[i,1],x[i,2]),cbind(rep(c("fq_bias","sq_bias","tq_bias"),times=n.ests),as.character(rep(1:n.ests,each=3))))
bias_res <- simres[,c("whichbeta","nsamp",bias_names),]
bias_fq_names <- sapply(1:n.ests,function(i,x) paste0("fq_bias",i))
bias_sq_names <- sapply(1:n.ests,function(i,x) paste0("sq_bias",i))
bias_tq_names <- sapply(1:n.ests,function(i,x) paste0("tq_bias",i))
bias_df <- reshape(bias_res,varying=list(bias_fq_names,bias_sq_names,bias_tq_names),direction="long")
names(bias_df)[3] <- "mod"
bias_df$bias <- character(length=nrow(bias_df))
bias_df$mod <- as.factor(bias_df$mod)
bias_df$mod_names <- "ssMOReTreeS Collapsed"
bias_df$mod_names[bias_df$mod==8] <- "ssMOReTreeS Individual"
bias_df$mod_names[bias_df$mod==2] <- "Uncollapsed"
bias_df$mod_names[bias_df$mod==3] <- "True Grouping"
bias_df$mod_names[bias_df$mod==4] <- "Adhoc Grouping 1"
bias_df$mod_names[bias_df$mod==5] <- "Adhoc Grouping 2"
bias_df$mod_names[bias_df$mod==6] <- "Adhoc Grouping 3"
bias_df$mod_names[bias_df$mod==7] <- "Fully Collapsed"
bias_df$whichbeta[bias_df$whichbeta=="Scenario 1"] <- "1"
bias_df$whichbeta[bias_df$whichbeta=="Scenario 2"] <- "2"
bias_df$mod_names <- factor(bias_df$mod_names,c("ssMOReTreeS Collapsed",
                                                    "ssMOReTreeS Individual",
                                                    "True Grouping",
                                                    "Uncollapsed",
                                                    "Adhoc Grouping 1",
                                                    "Adhoc Grouping 2",
                                                    "Adhoc Grouping 3",
                                                    "Fully Collapsed"))
digits <- 2
for(i in 1:nrow(bias_df)){
  frmt <- paste0("%.",digits,"f")
  bias_df$bias[i] <- paste0(sprintf(frmt,bias_df$sq_bias1[i])," (",sprintf(frmt,bias_df$fq_bias1[i]),",",sprintf(frmt,bias_df$tq_bias1[i]),")")
}
bias_df2 <- dcast(bias_df[,c("whichbeta","mod_names","nsamp","bias")],whichbeta + mod_names ~ nsamp, value.var="bias",fun.aggregate=)

require(xtable)
bias_xtable <- xtable(bias_df2)
write(print(bias_xtable,floating=FALSE,include.rownames = FALSE),file="./moretrees/tables_and_figures/tableA1.tex")

############ RMSE of group estimates for individual true betas ###########

n.ests <- 8
rmse_fq_names <- sapply(1:n.ests,function(i,x) paste0(x[i,1],x[i,2]),cbind(rep("fq_rmse",times=n.ests),as.character(1:n.ests)))
rmse_sq_names <- sapply(1:n.ests,function(i,x) paste0(x[i,1],x[i,2]),cbind(rep("sq_rmse",times=n.ests),as.character(1:n.ests)))
rmse_tq_names <- sapply(1:n.ests,function(i,x) paste0(x[i,1],x[i,2]),cbind(rep("tq_rmse",times=n.ests),as.character(1:n.ests)))
nclust_fq_names <- sapply(1:n.ests,function(i,x) paste0(x[i,1],x[i,2]),cbind(rep("fq_nclust",times=n.ests),as.character(1:n.ests)))
nclust_sq_names <- sapply(1:n.ests,function(i,x) paste0(x[i,1],x[i,2]),cbind(rep("sq_nclust",times=n.ests),as.character(1:n.ests)))
nclust_tq_names <- sapply(1:n.ests,function(i,x) paste0(x[i,1],x[i,2]),cbind(rep("tq_nclust",times=n.ests),as.character(1:n.ests)))
ARI_fq_names <- sapply(1:n.ests,function(i,x) paste0("fq_ARI",i))
ARI_sq_names <- sapply(1:n.ests,function(i,x) paste0("sq_ARI",i))
ARI_tq_names <- sapply(1:n.ests,function(i,x) paste0("tq_ARI",i))
all_names <- c(rmse_fq_names,rmse_sq_names,rmse_tq_names,nclust_fq_names,nclust_sq_names,nclust_tq_names,ARI_fq_names,ARI_sq_names,ARI_tq_names)
simres <- simres[,c(names(params),"nsamp_lab",all_names)]
simres.long <- reshape(simres,varying=list(rmse_fq_names,rmse_sq_names,rmse_tq_names,nclust_fq_names,nclust_sq_names,nclust_tq_names,ARI_fq_names,ARI_sq_names,ARI_tq_names),direction="long")
names(simres.long)[6] <- "mod"
simres.long$mod <- as.factor(simres.long$mod)
simres.long$mod_names <- "ssMOReTreeS Collapsed"
simres.long$mod_names[simres.long$mod==8] <- "ssMOReTreeS Individual"
simres.long$mod_names[simres.long$mod==2] <- "Uncollapsed"
simres.long$mod_names[simres.long$mod==3] <- "True Grouping"
simres.long$mod_names[simres.long$mod==4] <- "Adhoc Grouping 1"
simres.long$mod_names[simres.long$mod==5] <- "Adhoc Grouping 2"
simres.long$mod_names[simres.long$mod==6] <- "Adhoc Grouping 3"
simres.long$mod_names[simres.long$mod==7] <- "Fully Collapsed"
simres.long$mod_names <- factor(simres.long$mod_names,rev(c("ssMOReTreeS Collapsed",
                                                        "ssMOReTreeS Individual",
                                                        "True Grouping",
                                                        "Uncollapsed",
                                                        "Adhoc Grouping 1",
                                                        "Adhoc Grouping 2",
                                                        "Adhoc Grouping 3",
                                                        "Fully Collapsed")))

################# Figure 3- RMSE for different models #####################

rmse_plot <- ggplot(simres.long, aes(x=sq_rmse1, y=mod_names, col=mod)) + 
  geom_point() +
  geom_errorbarh(aes(xmin=fq_rmse1, xmax=tq_rmse1))+
  scale_x_continuous(trans="log10",breaks=c(0.0001,0.001,0.01,0.1,1,10,100,1000)) +
  facet_grid(nsamp_lab~whichbeta) +
  ylab("Model") +
  xlab("Median (Interquartile Range) RMSE (log 10 scale)") +
  theme_bw(base_size=16) +
  theme(legend.position="none")

pdf(file="./moretrees/tables_and_figures/figure3.pdf",width=9,height=8)
rmse_plot
dev.off()

################# Figure A2 - RMSE vs ARI and number of clusters #####################

plot_ARI_rmse <- ggplot(data=simres.long[simres.long$mod == 1,],aes(x=sq_ARI1,y=sq_rmse1,col=nsamp_lab)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=fq_rmse1, ymax=tq_rmse1, col=nsamp_lab), width=0.05,lwd=1) +
  geom_errorbarh(aes(xmin=fq_ARI1, xmax=tq_ARI1, col=nsamp_lab), height=0.002,lwd=1) +
  theme_bw(base_size=25) +
  geom_vline(xintercept=1,col="grey") +
  xlab("Median (Interquartile Range) ARI") +
  ylab("Median (Interquartile Range) RMSE") +
  guides(col=guide_legend(title="Sample Size")) +
  facet_wrap(whichbeta~.)

pdf(file="./moretrees/tables_and_figures/tableA2_a.pdf",width=12,height=6)
plot_ARI_rmse
dev.off()

# With nclust on x axis

plot_nclust_rmse <- ggplot(data=simres.long[simres.long$mod == 1,],aes(x=sq_nclust1,y=sq_rmse1,col=nsamp_lab)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=fq_rmse1, ymax=tq_rmse1, col=nsamp_lab), width=0.05,lwd=1) +
  geom_errorbarh(aes(xmin=fq_nclust1, xmax=tq_nclust1, col=nsamp_lab), height=0.002,lwd=1) +
  theme_bw(base_size=25) +
  geom_vline(xintercept=1,col="grey") +
  xlab("Median (Interquartile Range) Number of Clusters Estimated") +
  ylab("Median (Interquartile Range) RMSE") +
  guides(col=guide_legend(title="Sample Size")) +
  geom_vline(xintercept=simres.long$fq_nclust1[simres.long$mod_name == "True Grouping"][1],col="red") +
  facet_wrap(whichbeta~.)

pdf(file="./moretrees/tables_and_figures/tableA2_b.pdf",width=12,height=6)
plot_nclust_rmse
dev.off()
