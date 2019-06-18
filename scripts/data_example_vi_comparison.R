# ----------------------------------------------------------------- #
# ------------------- Data analysis: VI comparison ---------------- #
# ----------------------------------------------------------------- #

# direc <- "../moretrees/" # path of the moretrees repository
# direc <- "~/Documents/PhD_Papers/Air_pollution/R_code/MORETreeS/moretrees/"
direc <- "/nfs/home/E/ethomas/shared_space/ci3_nsaph/Emma/R_code/moretrees"
setwd(direc)

if(!dir.exists("./data_example_results")) dir.create("./data_example_results")

### Load functions ###
source("scripts/VI_functions.R")
source("scripts/processing_functions.R")
require(igraph)
require(Matrix)
require(doParallel)

### Load ICD9 tree ###
load("./simulation_inputs/inputs.Rdata")

######### Algorithm parameters #########

datArgs <- as.integer(as.character(commandArgs(trailingOnly = TRUE))) # Use to call arguments from the command line
# datArgs <- c(0,3,1E5,1E-8,10) # Alternatively, enter arguments directly in R

sim <- datArgs[1]+1 # which simulated dataset (integer from 1 to 10)
nrestarts <- datArgs[2] # number of random restarts for VI algorithm
m.max <- datArgs[3] # maximum number of time steps
tol <- datArgs[4] # tolerance for convergence
m.print <- datArgs[5] # how often to print progress for VI algorithm

######################## Load simulated data ########################

load(file=paste0("./simulation_inputs/simulate_data_comparison",sim,".Rdata"))
Y.sim <- sapply(Z.sim,length)

######################## Variational inference ########################

# Adhoc collapsing estimates
adhoc_coeffs <- adhoc_collapsing(Z.sim,Y.sim,pL,groups)

# Initial values for node coefficients
nodes_init <- initial_node_coeffs(Z.sim,Y.sim,uncollapsed=adhoc_coeffs[,1],p,pL,leaf.descendants,ancestors)

# For parallelization
registerDoParallel(cores=nrestarts)

# Run VI algorithm using nrestarts random restarts
restarts_vi <- foreach(j = 1:nrestarts) %dopar% {
  
  # saving output to track number of time steps
  pr <- paste0("./data_example_results/comparison_vi_print",sim,"_restart",j,".txt")
  sink(file=pr)
  
  # profiling to test speed of VI vs. MCMC
  prof <- paste0("./data_example_results/comparison_vi_prof",sim,"_restart",j,".out")
  Rprof(file=prof,memory.profiling=TRUE)
  
  # run VI algorithm
  out_vi <- VI_binary_ss(Z=Z.sim,Y=Y.sim,n=nsamp,p=p,pL=pL,ancestors=ancestors,
                         leaf.descendants=leaf.descendants,cutoff=0.5,mu_gamma_init=nodes_init,
                         tol=tol,m.max=m.max,m.print=m.print,more=FALSE,update_hyper=T,update_hyper_freq=10)
  
  Rprof(NULL) # close Rprof
  print(summaryRprof(prof,lines="hide",memory="both")) # summarize Rprof results
  file.remove(prof) # need to remove profile files as they are very large
  
  sink() # close sink
  
  out_vi
  
}

# Choose final model with highest ELBO
ELBOS <- numeric(nrestarts)
for(j in 1:nrestarts){
  ELBOS[j] <- restarts_vi[[j]]$ELBO
}
out_vi <- restarts_vi[[which.max(ELBOS)]]

# Save VI results
res <- paste0("./data_example_results/comparison_vi_results",sim,".Rdata")
save(out_vi,file=res)
