# ------------------------------------------------------- #
# ------- Data analysis - cross validation -------------- #
# ------------------------------------------------------- #

# direc <- "../moretrees/" # path of the moretrees repository
direct <- "/Users/emt380/Documents/PhD_Papers/Air_pollution/R_code/MORETreeS/moretrees/"
setwd(direc)

#### Create directory for saving results ###
if(!dir.exists("./data_example_results")) dir.create("./data_example_results")

##### Load functions ######
source("./scripts/VI_functions.R")
source("./scripts/processing_functions.R")
require(igraph)
require(Matrix)

### load ICD9 tree ###
load("./simulation_inputs/inputs.Rdata")
# Extract list of relevant ICD9 codes
codes <- names(V(tree)[V(tree)$leaf])

######### Algorithm parameters #########

datArgs <- as.integer(as.character(commandArgs(trailingOnly = TRUE))) # Use to call arguments from the command line
# datArgs <- c(0,1E5,1E-8) # Alternatively, enter arguments directly in R

fold <- datArgs[1]+1 # which fold for cv (integer from 1 to 10)
m.max <- datArgs[2] # maximum number of time steps
tol <- datArgs[3] # tolerance for convergence

# UNCOMMENT THE BELOW BEFORE MERGE WITH MAIN BRANCH

# ############### Prepare data ###############
# 
# # Load data
# load(file="data/moretrees_CC_data.Rdata")
# # Load cv folds
# load(file="data/cv_folds.Rdata")
# 
# ############### Out-of-sample prediction via 10-fold CV ###############
# 
# # Function for computing log-likelihood component for each outcome
# ll.fun <- function(v,beta,Z.test){ 
#   sum(loglogit(beta[v]*Z.test[[v]]))
# } 
# 
# # Extract training and test data
# Y.train <- numeric(length=pL)
# Z.train <- list()
# Y.test <- numeric(length=pL)
# Z.test <- list()
# for(v in 1:pL){
#   Y.test[v] <- sum(folds[[v]] == fold)
#   Y.train[v] <- Y[v] - Y.test[v]
#   Z.test[[v]] <- Z[[v]][folds[[v]]==fold]
#   Z.train[[v]] <- Z[[v]][folds[[v]]!=fold]
# }
# # Discard original data to save on memory
# rm(Y,Z)
# # Adhoc collapsing estimates
# adhoc_coeffs_fold <- adhoc_collapsing(Z.train,Y.train,pL,groups)
# # Initial values for node coefficients
# nodes_init_fold <- initial_node_coeffs(Z.train,Y.train,uncollapsed=adhoc_coeffs_fold[,1],p,pL,leaf.descendants,ancestors)
# # Run ssMOReTreeS
# mod <- VI_binary_ss(Z=Z.train,Y=Y.train,n=sum(Y.train),p=p,pL=pL,ancestors=ancestors,
#                     leaf.descendants=leaf.descendants,cutoff=0.5,mu_gamma_init=nodes_init_fold,
#                     tol=tol,m.max=m.max,m.print=10,more=FALSE,update_hyper=TRUE,update_hyper_freq=10)
# # Get test set log likelihoods
# ll.moretrees.group <- sum(sapply(1:pL,ll.fun,beta=mod$moretrees_est,Z=Z.test))/sum(Y.test)
# beta.indiv <- indiv.beta.calc(VI_params=mod$VI_params,ancestors=ancestors,pL=pL,p=p)
# ll.moretrees.indiv <- sum(sapply(1:pL,ll.fun,beta=beta.indiv,Z=Z.test))/sum(Y.test)
# ll.adhoc <- numeric(ncol(adhoc_coeffs_fold))
# for(i in 1:ncol(adhoc_coeffs_fold)){
#   ll.adhoc[i] <- sum(sapply(1:pL,ll.fun,beta=adhoc_coeffs_fold[,i],Z=Z.test))/sum(Y.test)
# }
# # Result
# ll.cv <- as.data.frame(matrix(c(fold,ll.moretrees.group,ll.moretrees.indiv,ll.adhoc),nrow=1))
# names(ll.cv) <- c("fold","moretrees_groups","moretrees_indiv","uncollapsed","sim_groups","adhoc1","adhoc2","adhoc3","fully_collapsed")

############### Testing values of tuning parameter ###############
# loading results of VI algo- DELETE THIS LATER
load(paste0("./data_example_results/data_example_cv_fold",fold,".Rdata"))

# Values of the tuning parameter to test
tp <- seq(0,1,0.01)

# Get ancestor matrix A
A <- t(as_adj(tree,sparse = T))
A <- expm(A)
A[A>0] <- 1 

# Extract relevant VI params
mu_gamma <- mod$mu_gamma
p_vi <- exp(loglogit(mod$VI_params$u_s))

ll.tp <- numeric(length(tp))
for(i in length(tp_vals)){
  beta.t <- beta <- A %*% (mu_gamma * (p_vi >= t[i]))
  ll.tp[i] <- sum(sapply(1:pL,ll.fun,beta=beta.t,Z=Z.test))/sum(Y.test)
}


############### Save results ###############

save(ll.cv,mod,ll.tp,tp,file=paste0(direc,"data_example_results/data_example_cv_fold",fold,".Rdata"))
