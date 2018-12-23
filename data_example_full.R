# ------------------------------------------------------- #
# ------------------ Data analysis ---------------------- #
# ------------------------------------------------------- #

direc <- "/nfs/home/E/ethomas/shared_space/ci3_nsaph/Emma/R_code/MORETreeS/moretrees/"
#direc <- "/Users/emt380/Documents/PhD_Papers/Air_pollution/R_code/MORETreeS/moretrees/"
setwd(direc)

if(!dir.exists("./data_example_results")) dir.create("./data_example_results")

##### Load functions ######
source("VI_functions.R")
source("processing_functions.R")
require(igraph)

### load ICD9 tree ###
load("./simulation_inputs/inputs.Rdata")
# Extract list of relevant ICD9 codes
codes <- names(V(tree)[V(tree)$leaf])

######### Algorithm parameters #########

datArgs <- as.integer(as.character(commandArgs(trailingOnly = TRUE)))
#datArgs <- c(3,1E5,1E-8,10)

nrestarts <- datArgs[1] # number of random restarts
m.max <- datArgs[2] # maximum number of time steps
tol <- datArgs[3] # tolerance for convergence

############### Prepare data ###############

# Load data
load(file="/nfs/home/E/ethomas/shared_space/ci3_nsaph/Emma/Data/Case_crossover_data/moretrees_CC_data.Rdata")

# Set seed
set.seed(4958764)

############### Run analysis on wole dataset ###############

source("./Master/functions_MORETreeS_VI_case_control2.R")

# Adhoc collapsing estimates
adhoc_coeffs <- adhoc_collapsing(Z,Y,pL,groups)

# Initial values for node coefficients
nodes_init <- initial_node_coeffs(Z,Y,uncollapsed=adhoc_coeffs[,1],p,pL,leaf.descendants,ancestors)

# For parallelization
require(doParallel)
registerDoParallel(cores=nrestarts)

# Run spike & slab model using nrestarts random restarts
restarts_ss <- foreach(j = 1:nrestarts) %dopar% {
  
  # out_ss <- 
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

# Run usual ML logistic regression model using groups discovered by moretrees
beta_est <- final_ss$moretrees_est
groups <- as.numeric(as.factor(beta_est))
beta.ml.groups <- data.frame(beta_est=numeric(max(groups)),beta_cil=numeric(max(groups)),beta_ciu=numeric(max(groups)))
for(g in 1:max(groups)){
  which.dat <- groups==g
  Y.g <- sum(Y[which.dat])
  if(Y.g == 0){
    beta.groups[which.dat,i] <- NA
  } else {
    beta_ml <- glm(rep(1,Y.g) ~ 0 + unlist(Z[which.dat]),family="binomial")
    beta.ml.groups[g,] <- c(beta_ml$coefficients[1],
                          beta_ml$coefficients[1]+qnorm(0.025)*sqrt(vcov(beta_ml)),
                          beta_ml$coefficients[1]+qnorm(0.975)*sqrt(vcov(beta_ml)))
  }
}

############### Save results ###############

save(beta_est,groups,beta.ml.groups,Y,final_ss,adhoc_coeffs,ELBOS,file = paste0(direc,"data_example_results/data_example_full.Rdata"))
