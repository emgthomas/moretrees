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
source("scripts/logit_spike_edit.R")
require(igraph)
require(Matrix)
require(BoomSpikeSlab)

######### Algorithm parameters #########

datArgs <- as.integer(as.character(commandArgs(trailingOnly = TRUE))) # Use to call arguments from the command line
# datArgs <- c(1,3,10,1,1) # Alternatively, enter arguments directly in R
runagain <- c(1,9) # Delete later; rerunning some chains that didn't work
datArgs[1] <- runagain[datArgs[1]+1]

nchains <- datArgs[2] # how may parallel chains to run
sim <- datArgs[1] %/% nchains + 1 # which simulated dataset (integer from 1 to 10)
chain <- datArgs[1] %% nchains + 1 # which chain (integer from 1 to nchains)
niter <- datArgs[3] # number of MCMC samples
nthreads <- datArgs[4] # number of threads for data augmentation (see ?logit.spike)
ping <- datArgs[5] # print progress report after ping samples

######################## Load data ########################

# Load ICD9 tree
load("./simulation_inputs/inputs.Rdata")

# Load VI results
VI_res <- paste0("./data_example_results/comparison_vi_results",sim,".Rdata")
load(VI_res)

# Load simulated data
load(file=paste0("./simulation_inputs/simulate_data_comparison",1,".Rdata"))

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
Xmat <- bdiag(Z.sim)
Xmat <- cbind(Matrix(0,nrow=nrow(Xmat),ncol=p-pL,sparse=T),Xmat)
Xstar <- Xmat %*% A

# dummy outcome
Yvec <- rep(1,nrow(Xstar))

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

# saving output to track number of time steps
pr <- paste0("./data_example_results/comparison_mcmc_print",sim,"_chain",chain,".txt")
sink(file=pr)

# checking if we already have some iterations to work with
res <- paste0("./data_example_results/comparison_mcmc_samples",sim,"_chain",chain,".csv")
if(exists(res)){
  nL <- countLines(res)
  gamma0 <- as.numeric(read.csv(res, header=FALSE, skip=nL-1)) # Read only last line to restart from most recent iteration
} else {
  # If not, use VI results as initial values
  gamma0 <- gamma_vi
}

# profiling to test speed of VI vs. MCMC
prof <- paste0("./data_example_results/comparison_mcmc_prof",sim,"_chain",chain,".out")
# Use append=TRUE in case we have already run some iterations
Rprof(file=prof,memory.profiling=TRUE,append=TRUE)

# run MCMC
out_mcmc <- logit.spike.edit(Y ~ 0 + ., data=data.frame(Y=Yvec,X=as.matrix(Xstar)),
                             niter=niter,
                             prior=prior,
                             initial.value=gamma0,
                             nthreads=nthreads,
                             ping=ping)

Rprof(NULL) # close Rprof
print(summaryRprof(prof,lines="hide",memory="both")) # summarize Rprof results
file.remove(prof) # need to remove profile files as they are very large

######################## Save MCMC samples #######################
sgamma_mcmc <- out_mcmc$beta
colnames(sgamma_mcmc) <- colnames(Xstar)

write.table(sgamma_mcmc,file=res,sep=",",col.names=F,row.names=F,append=TRUE)

sink() # close sink
