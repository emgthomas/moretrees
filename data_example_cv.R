# ------------------------------------------------------- #
# ------- Data analysis - cross validation -------------- #
# ------------------------------------------------------- #

direc <- "/nfs/home/E/ethomas/shared_space/ci3_nsaph/Emma/R_code/MORETreeS/"
# direc <- "/Users/emt380/Documents/PhD_Papers/Air_pollution/R_code/MORETreeS/"
setwd(direc)

##### Load functions ######
# source("./Master/functions_MORETreeS_VI_case_control2.R")
source("./moretrees/VI_functions.R")
require(igraph)
indiv.beta.calc <- function(VI_params,ancestors,pL,p){
  pi_s <- 1/(1+exp(-VI_params$u_s))
  mu_gamma <- VI_params$mu_gamma
  beta_est <- numeric(length=pL)
  for(v in 1:pL){
    u <- ancestors[[v+p-pL]]
    beta_est[v] <- sum(mu_gamma[u]*pi_s[u])
  }
  return(beta_est)
}

### load ICD9 tree ###
load("./Data/sim_trees_CC.Rdata")
tree <- tree[[2]]
p <- p[2]
pL <- pL[2]
ancestors <- ancestors[[2]]
leaf.descendants <- leaf.descendants[[2]]
n.leaf.descendants <- n.leaf.descendants[2]
groups <- groups[[2]]
# Extract list of relevant ICD9 codes
codes <- names(V(tree)[V(tree)$leaf])

######### Algorithm parameters #########

datArgs <- as.integer(as.character(commandArgs(trailingOnly = TRUE)))
#datArgs <- c(10,1E5,1E-8)

fold <- datArgs[1]+1 # which fold for cv
m.max <- datArgs[2] # maximum number of time steps
tol <- datArgs[3] # tolerance for convergence

############### Prepare data ###############

# Load data
load(file="/nfs/home/E/ethomas/shared_space/ci3_nsaph/Emma/Data/Case_crossover_data/moretrees_CC_data.Rdata")

# Load cv folds
load(file="./Data/cv_folds.Rdata")

source("./Master/functions_MORETreeS_VI_case_control2.R")

############### GoF via 10-fold CV ###############

# Function for computing log-likelihood component for each outcome
ll.fun <- function(v,beta,Z.test){ 
  sum(loglogit(beta[v]*Z.test[[v]]))
} 

# Extract training and test data
Y.train <- numeric(length=pL)
Z.train <- list()
Y.test <- numeric(length=pL)
Z.test <- list()
for(v in 1:pL){
  Y.test[v] <- sum(folds[[v]] == fold)
  Y.train[v] <- Y[v] - Y.test[v]
  Z.test[[v]] <- Z[[v]][folds[[v]]==fold]
  Z.train[[v]] <- Z[[v]][folds[[v]]!=fold]
}
# Discard original data to save on memory
rm(Y,Z)
# Adhoc collapsing estimates
adhoc_coeffs_fold <- adhoc_collapsing(Z.train,Y.train,pL,groups)
# Initial values for node coefficients
nodes_init_fold <- initial_node_coeffs(Z.train,Y.train,uncollapsed=adhoc_coeffs_fold[,1],p,pL,leaf.descendants,ancestors)
# Run ssMOReTreeS
mod <- VI_binary_ss(Z=Z.train,Y=Y.train,n=sum(Y.train),p=p,pL=pL,ancestors=ancestors,
                    leaf.descendants=leaf.descendants,cutoff=0.5,mu_gamma_init=nodes_init_fold,
                    tol=tol,m.max=m.max,m.print=10,more=FALSE,update_hyper=TRUE,update_hyper_freq=10)
# Get test set log likelihoods
ll.moretrees.group <- sum(sapply(1:pL,ll.fun,beta=mod$moretrees_est,Z=Z.test))/sum(Y.test)
beta.indiv <- indiv.beta.calc(VI_params=mod$VI_params,ancestors=ancestors,pL=pL,p=p)
ll.moretrees.indiv <- sum(sapply(1:pL,ll.fun,beta=beta.indiv,Z=Z.test))/sum(Y.test)
ll.adhoc <- numeric(ncol(adhoc_coeffs_fold))
for(i in 1:ncol(adhoc_coeffs_fold)){
  ll.adhoc[i] <- sum(sapply(1:pL,ll.fun,beta=adhoc_coeffs_fold[,i],Z=Z.test))/sum(Y.test)
}
# Result
ll.cv <- as.data.frame(matrix(c(fold,ll.moretrees.group,ll.moretrees.indiv,ll.adhoc),nrow=1))
names(ll.cv) <- c("fold","moretrees_groups","moretrees_indiv","uncollapsed","sim_groups","adhoc1","adhoc2","adhoc3","fully_collapsed")

############### Save results ###############

save(ll.cv,mod,file=paste0(direc,"Results/Data_Example/data_example2_cv_fold",fold,".Rdata"))
#write.table(ll.cv,file=paste0(direc,"Results/Data_Example/data_example_cv_fold",fold,".csv"))
