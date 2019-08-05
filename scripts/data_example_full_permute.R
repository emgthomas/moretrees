# --------------------------------------------------------------------------------- #
# ------------------ Data analysis: permuting leaves of tree ---------------------- #
# --------------------------------------------------------------------------------- #

# direc <- "../moretrees/" # path of the moretrees repository
# direc <- "/Users/emt380/Documents/PhD_Papers/Air_pollution/R_code/MORETreeS/moretrees/"
direc <- "/nfs/home/E/ethomas/shared_space/ci3_nsaph/Emma/R_code/moretrees"
setwd(direc)

if(!dir.exists("./data_example_results")) dir.create("./data_example_results")

##### Load functions ######
source("scripts/VI_functions.R")
source("scripts/processing_functions.R")
require(igraph)

### load ICD9 tree ###
load("./simulation_inputs/inputs.Rdata")
# Extract list of relevant ICD9 codes
codes <- names(V(tree)[V(tree)$leaf])

######### Input parameters #########

datArgs <- as.integer(as.character(commandArgs(trailingOnly = TRUE))) # Use to call arguments from the command line
# datArgs <- c(0,3,1E5,1E-16) # Alternatively, enter arguments directly in R

perm <- datArgs[1]+1 # permuation number
nrestarts <- datArgs[2] # number of random restarts
m.max <- datArgs[3] # maximum number of time steps
tol <- datArgs[4] # tolerance for convergence

############### Prepare data ###############

# Load data
load(file="/nfs/home/E/ethomas/shared_space/ci3_nsaph/Emma/Data/moretrees_data/moretrees_CC_data.Rdata")

# Load permutations
load(file="./simulation_inputs/permutations.Rdata")

# Load original results
load(file="./data_example_results/data_example_full.Rdata")

# Permute data
set.seed(permutations[perm])
for(v in 1:pL){
  Z[[v]] <- Z[[v]]*sample(c(-1,1),Y[v],replace=T)
}

############### Run analysis on whole dataset ###############

# Adhoc collapsing estimates
adhoc_coeffs <- adhoc_collapsing(Z,Y,pL,groups)

# Initial values for node coefficients
nodes_init <- initial_node_coeffs(Z,Y,uncollapsed=adhoc_coeffs[,1],p,pL,leaf.descendants,ancestors)

# For parallelization
require(doParallel)
registerDoParallel(cores=nrestarts)

# Run spike & slab model using nrestarts random restarts
restarts_ss <- foreach(j = 1:nrestarts) %dopar% {
  
  #out_ss <- 
  VI_binary_ss(Z=Z,Y=Y,n=sum(Y),p=p,pL=pL,ancestors=ancestors,
               leaf.descendants=leaf.descendants,cutoff=0.5,mu_gamma_init=nodes_init,
               tol=tol,m.max=m.max,m.print=10,more=FALSE,update_hyper=TRUE,update_hyper_freq=10)
  
}

# Choose final model with highest ELBO
ELBOS <- numeric(nrestarts)
for(j in 1:nrestarts){
  ELBOS[j] <- restarts_ss[[j]]$ELBO
}
final_ss <- restarts_ss[[which.max(ELBOS)]]

# Get final estimates
beta_est <- final_ss$moretrees_est
groups <- as.numeric(as.factor(beta_est))

############### Save results ###############

save(beta_est,groups,final_ss,ELBOS,
     file = paste0("./data_example_results/data_example_full_perm",perm,".Rdata"))
