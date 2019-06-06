# -------------------------------------------------------------------- #
# --------------------------- VI for MORETREES ----------------------- #
# ------------- Case-control setting with spike & slab prior --------- #
# -------------------------------------------------------------------- #

direc <- "~/Documents/PhD_Papers/Air_pollution/R_code/MORETreeS/moretrees/" # set to path of the moretrees repository
setwd(direc)

##### Load functions ######
source("scripts/VI_functions.R")
require(igraph)

######################## SET UP TREE ########################

### load test tree ###
load("simulation_inputs/inputs.Rdata")

######### Simulation parameters #########

# simArgs <- as.integer(as.character(commandArgs(trailingOnly = TRUE))) # Use to call arguments from the command line
simArgs <- c(0,10,1,1E5,1E-8,2) # Alternatively, enter arguments directly in R

simidx <- simArgs[1] + 1 # which set of parameters to use. Can take integer values 1 through 8.
nsims <- simArgs[2] # number of sims
nrestarts <- simArgs[3] # number of random restarts
m.max <- simArgs[4] # maximum number of time steps
tol <- simArgs[5] # tolerance for convergence
ncores <- simArgs[6] # number of cores

# Retrieve simulation parameters using process number
nsamp <- params[simidx,]$nsamp # sample size
whichbeta <- params[simidx,]$whichbeta # which set of betas to use: 1 (large beta) or 2 (smaller betas based on real data)

# Retrieve inputs based on simulation parameters
if(whichbeta==1){
  beta <- V(tree)$beta1[(p-pL+1):p]
} else {
  beta <- V(tree)$beta2[(p-pL+1):p]
}

######################## PREPARE FILES ########################

# For parallelization
require(doParallel)
registerDoParallel(cores=ncores)

# For storing results
path_file <- "simulation_results/"
sim_params <- paste("n",nsamp,"_beta",whichbeta,sep="")
params$nsamp2 <- as.character(params$nsamp)
params$nsamp2[params$nsamp==1E05] <- "1e\\+05"
params$nsamp2[params$nsamp==1E06] <- "1e\\+06"
nsamp2 <- params$nsamp2[simidx]
sim_params2 <- paste("n",nsamp2,"_beta",whichbeta,sep="")

# Clear existing partial result files
beta_files <- list.files(path=path_file,pattern=paste("part_betasims_",sim_params2,"*",sep=""))
sapply(beta_files, function(fn) if (file.exists(fn)) file.remove(fn))
VI_files <- list.files(path=path_file,pattern=paste("part_VIsims_",sim_params2,"*",sep=""))
sapply(VI_files, function(fn) if (file.exists(fn)) file.remove(fn))
hyper_files <- list.files(path=path_file,pattern=paste("part_hyperparams_",sim_params2,"*",sep=""))
sapply(hyper_files, function(fn) if (file.exists(fn)) file.remove(fn))

# Creating log file
logfile <- paste0(path_file,"log_n",nsamp,"_beta",whichbeta,".csv")
file.create(logfile)

######################## LOAD DATA ########################

# Load Z.sim and Y.sim- this depends on sample size
load(file=paste0("simulation_inputs/simdat_n",nsamp2,".Rdata"))

######################## RUN SIMS ########################

print("Running simulations...")

foreach(i=1:nsims,.combine=rbind,.errorhandling="remove") %dopar% {
  
  # Simulate case/control status
  for(v in 1:pL){
    if(Y.sim[v] > 0){
      beta.v <- beta[v]
      # Compute probability unit 1 is a case
      prob_unit1 <- 1/(1+exp(-beta.v*Z.sim[[v]]))
      probs <- as.data.frame(rbind(prob_unit1,1-prob_unit1))
      # Randomly assign case or control status
      cc <- sapply(probs,sample,x=c(1,-1),size=1,replace=FALSE)
      Z.sim[[v]] <- Z.sim[[v]]*cc
    }
  }
  
  # Adhoc collapsing estimates
  adhoc_coeffs <- adhoc_collapsing(Z.sim,Y.sim,pL,groups)
  
  # Initial values for node coefficients
  nodes_init <- initial_node_coeffs(Z.sim,Y.sim,uncollapsed=adhoc_coeffs[,1],p,pL,leaf.descendants,ancestors)
  
  # Run ssMOReTreeS model
  out_ss <- VI_binary_ss(Z=Z.sim,Y=Y.sim,n=nsamp,p=p,pL=pL,ancestors=ancestors,
                         leaf.descendants=leaf.descendants,cutoff=0.5,mu_gamma_init=nodes_init,
                         tol=tol,m.max=m.max,m.print=m.max+1,more=FALSE,update_hyper=T,update_hyper_freq=10)
  
  write.table(paste("Simulation",i,"completed",sep=" "), file = logfile, row.names =FALSE, col.names = FALSE,sep = ",", append = TRUE)
  
  outfile_beta <- file(paste(path_file,"part_betasims_",sim_params,"_",Sys.info()[['nodename']], Sys.getpid(),".csv",sep=""), open="a")
  write.table(cbind(rep(i,pL),out_ss$moretrees_est,adhoc_coeffs), file = outfile_beta, row.names =FALSE, col.names = FALSE,sep = ",", append = TRUE)
  close(outfile_beta)
  
  outfile_VI <- file(paste(path_file,"part_VIsims_",sim_params,"_",Sys.info()[['nodename']], Sys.getpid(),".csv",sep=""), open="a")
  write.table(cbind(rep(i,p),out_ss$VI_params), file = outfile_VI, row.names =FALSE, col.names = FALSE,sep = ",", append = TRUE)
  close(outfile_VI)
  
  outfile_hyperparams <- file(paste(path_file,"part_hyperparams_",sim_params,"_",Sys.info()[['nodename']], Sys.getpid(),".csv",sep=""), open="a")
  write.table(rbind(c(i,out_ss$hyperparams)), file = outfile_hyperparams, row.names =FALSE, col.names = FALSE,sep = ",", append = TRUE)
  close(outfile_hyperparams)
  
  outfile_reached_max <- file(paste(path_file,"part_reached_max_",sim_params,"_",Sys.info()[['nodename']], Sys.getpid(),".csv",sep=""), open="a")
  write.table(rbind(c(i,out_ss$reached.max)), file = outfile_reached_max, row.names =FALSE, col.names = FALSE,sep = ",", append = TRUE)
  close(outfile_reached_max)
}

### Append results together and remove temporary files ###

# betas
beta_files <- list.files(path=path_file,pattern=paste("part_betasims_",sim_params2,"_*",sep=""))
outfile_beta <- paste0(path_file,"betasims_",sim_params,".csv")
file.create(outfile_beta)
write.table(rbind(c("sim","moretrees_est","uncollapse","truth","adhoc1","adhoc2","adhoc3","adhoc4")), file = outfile_beta, row.names =FALSE, col.names = FALSE,sep = ",", append = TRUE)
for(fn in beta_files){
  sims <- read.csv(file=paste0(path_file,fn),header=F,row.names=NULL)
  write.table(sims, file=outfile_beta, row.names =FALSE, col.names = FALSE,sep = ",", append = TRUE)
  file.remove(paste0(path_file,fn))
}

# VI params
VI_files <- list.files(path=path_file,pattern=paste("part_VIsims_",sim_params2,"_*",sep=""))
outfile_VI <- paste0(path_file,"VIsims_",sim_params,".csv")
file.create(outfile_VI)
write.table(rbind(c("sim","mu_gamma","sigma2_gamma","u_s")), file = outfile_VI, row.names=FALSE, col.names=FALSE, sep=",",append=TRUE)
for(fn in VI_files){
  sims <- read.csv(file=paste0(path_file,fn),header=F,row.names=NULL)
  write.table(sims, file=outfile_VI, row.names =FALSE, col.names = FALSE,sep = ",", append = TRUE)
  file.remove(paste0(path_file,fn))
}

# hyperparams
hyper_files <- list.files(path=path_file,pattern=paste("part_hyperparams_",sim_params2,"_*",sep=""))
outfile_hyper <- paste0(path_file,"hyperparams_",sim_params,".csv")
file.create(outfile_hyper)
write.table(rbind(c("sim","rho","tau")), file = outfile_hyper, row.names=FALSE, col.names=FALSE, sep=",",append=TRUE)
for(fn in hyper_files){
  sims <- read.csv(file=paste0(path_file,fn),header=F,row.names=NULL)
  write.table(sims, file=outfile_hyper, row.names =FALSE, col.names = FALSE,sep = ",", append = TRUE)
  file.remove(paste0(path_file,fn))
}

# reached max number of iterations
reached_max_files <- list.files(path=path_file,pattern=paste("part_reached_max_",sim_params2,"_*",sep=""))
outfile_reached_max <- paste0(path_file,"reached_max_",sim_params,".csv")
file.create(outfile_reached_max)
write.table(rbind(c("sim","reached_max")), file = outfile_reached_max, row.names=FALSE, col.names=FALSE, sep=",",append=TRUE)
for(fn in reached_max_files){
  sims <- read.csv(file=paste0(path_file,fn),header=F,row.names=NULL)
  write.table(sims, file=outfile_reached_max, row.names =FALSE, col.names = FALSE,sep = ",", append = TRUE)
  file.remove(paste0(path_file,fn))
}

# Done.

write.table("Done.", file = logfile, row.names =FALSE, col.names = FALSE,sep = ",", append = TRUE)

print("Done.")
