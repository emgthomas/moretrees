# ----------------------------------------------------------------- #
# ------------------ Data analysis: MCMC comparison --------------- #
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
require(BoomSpikeSlab)
require(doParallel)

### Load ICD9 tree ###
load("./simulation_inputs/inputs.Rdata")

### Load VI results ###

VI_res <- paste0("./data_example_results/comparison_vi_results",fold,".Rdata")
load(VI_res)

######### Algorithm parameters #########

datArgs <- as.integer(as.character(commandArgs(trailingOnly = TRUE))) # Use to call arguments from the command line
# datArgs <- c(0,3,10,1,1) # Alternatively, enter arguments directly in R

fold <- datArgs[1]+1 # which fold of data (integer from 1 to 10)
nchains <- datArgs[2] # how may parallel chains to run
niter <- datArgs[3] # number of MCMC samples
nthreads <- datArgs[4] # number of threads for data augmentation (see ?logit.spike)
ping <- datArgs[5] # print progress report after ping samples

######################## Load data ########################

# Load data
load(file="/nfs/home/E/ethomas/shared_space/ci3_nsaph/Emma/Data/moretrees_data/moretrees_CC_data.Rdata")
# Load folds
load(file="/nfs/home/E/ethomas/shared_space/ci3_nsaph/Emma/Data/moretrees_data/cv_folds.Rdata")

# Extract training data
for(v in 1:pL){
  Y[v] <- sum(folds[[v]] == fold)
  Z[[v]] <- Z[[v]][folds[[v]]==fold]
}

######################## Prepare data for MCMC ########################

# Get ancestor matrix A
A <- t(as_adj(tree,sparse = T))
A <- expm(A)
A[A>0] <- 1 

cat("\n\nCheck A is correct ancestor matrix\n\n")
A_check <- numeric(p)
for(v in 1:p) A_check[v] <- setequal(which(A[v,]>0),ancestors[[v]])
sum(A_check) == p

# Construct design matrix
Xmat <- bdiag(Z)
Xmat <- cbind(Matrix(0,nrow=nrow(Xmat),ncol=p-pL,sparse=T),Xmat)
Xstar <- Xmat %*% A
# # Do above matrix multiplication in chunks due to large matrix size
# Xstar <- Matrix(data=0,nrow=0,ncol=ncol(A),sparse=T)
# n <- nrow(Xmat)
# nchunks <- 20
# chunk_size <- round(n/nchunks)
# for(i in 1:nchunks){
#   if(i < nchunks){
#     idx <- ((i-1)*chunk_size+1):(i*chunk_size)
#   } else {
#     idx <- ((i-1)*chunk_size+1):n
#   }
#   Xstar <- rbind(Xstar,Xmat[idx,] %*% A)
#   print(i)
# }

# dummy outcome
Yvec <- rep(1,sum(Y))

# hyperparameters
rho <- out_vi$hyperparams[1]
tau <- out_vi$hyperparams[2]

# intial values for gamma_v*s_v
p_vi <- exp(loglogit(out_vi$VI_params$u_s))
s_vi <- p_vi >=0.5
gamma_vi <- out_vi$VI_params$mu_gamma*s_vi

# set up prior
prior <- LogitZellnerPrior(predictors=diag(rep(sqrt(p/tau)*2,p)),
                           prior.success.probability = 0.5,
                           expected.model.size=rho*p,
                           prior.information.weight = 1,
                           diagonal.shrinkage = 1,
                           max.flips=-1,
                           prior.inclusion.probabilities=rep(rho,p))

cat("\n\nCheck prior looks correct\n\n")
isDiagonal(prior$siginv)
sum(diag(prior$siginv)==1/tau) == p
sum(prior$prior.inclusion.probabilities == rho) == p

####################### Run MCMC #######################

# For parallelization
registerDoParallel(cores=nchains)

# Run spike & slab model using nchains chains
samples_mcmc <- foreach(j = 1:nchains) %dopar% {
  
  # saving output to track number of time steps
  pr <- paste0("./data_example_results/comparison_mcmc_print",fold,"_chain",j,".txt")
  sink(file=pr)
  
  # profiling to test speed of VI vs. MCMC
  prof <- paste0("./data_example_results/comparison_mcmc_prof",fold,"_chain",j,".out")
  Rprof(file=prof,memory.profiling=TRUE)
  
  # run MCMC
  out_mcmc <- logit.spike(Y ~ 0 + ., data=data.frame(Y=Yvec,X=as.matrix(Xstar)),
                          niter=niter,
                          prior=prior,
                          initial.value=gamma_vi,
                          nthreads=nthreads,
                          ping=ping)
  
  Rprof(NULL) # close Rprof
  print(summaryRprof(prof,lines="hide",memory="both")) # summarize Rprof results
  
  ######################## Save MCMC samples #######################
  sgamma_mcmc <- out_mcmc$beta
  colnames(sgamma_mcmc) <- colnames(Xstar)
  
  res <- paste0("./data_example_results/comparison_mcmc_samples",fold,"_chain",j,".csv")
  write.csv(sgamma_mcmc,file=res,row.names=F)
  
  sink() # close sink
}


