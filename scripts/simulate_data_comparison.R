# -------------------------------------------------------------------------- #
# ------------------- Simulated data for MCMC/VI comparison ---------------- #
# -------------------------------------------------------------------------- #

# direc <- "../moretrees/" # path of the moretrees repository
# direc <- "~/Documents/PhD_Papers/Air_pollution/R_code/MORETreeS/moretrees/"
direc <- "/nfs/home/E/ethomas/shared_space/ci3_nsaph/Emma/R_code/moretrees"
setwd(direc)

### Load functions and packages
source("scripts/processing_functions.R")

### Load ICD9 tree
load("./simulation_inputs/inputs.Rdata")

### Load data example results
load("./data_example_results/data_example_full.Rdata")

### Load data
load(file="/nfs/home/E/ethomas/shared_space/ci3_nsaph/Emma/Data/moretrees_data/moretrees_CC_data.Rdata")

### Sample size
n <- 1E6 # Size of simulated dataset
n_data <- sum(Y) # Size of original dataset
m <- sqrt(n_data/n) # Scaling factor for coefficients

### Compute individual beta estimates to use for simulation
beta_indiv_est <- indiv.beta.ci.calc(final_ss$VI_params,ancestors,pL,p)
# Scale up the beta estimates to account for smaller sample size
beta_sim <- beta_indiv_est$beta_indiv*m

### Simulate ten datasets ###
set.seed(48549)
for(i in 1:10){  
  
  # Simulate number of cases for each disease
  Y.sim <- 1:pL # ensure we have at least one case for each disease
  Y.sim <- c(Y.sim,sample(1:pL,size=n-pL,prob=Y/n_data,replace=T)) # sample the rest from true data
  Y.sim <- as.integer(table(sort(Y.sim)))
  
  # Simulate exposures and case/control status
  Z.sim <- list()
  for(v in 1:pL){
    # Sample exposures
    Z.sim[[v]] <- sample(Z[[v]],Y.sim[v],replace=T)
    # Compute probability unit 1 is a case
    prob_unit1 <- 1/(1+exp(-beta_sim[v]*Z.sim[[v]]))
    probs <- as.data.frame(rbind(prob_unit1,1-prob_unit1))
    # Randomly assign case or control status
    cc <- sapply(probs,sample,x=c(1,-1),size=1,replace=FALSE)
    Z.sim[[v]] <- Z.sim[[v]]*cc
  }
  
  # Save result
  save(Z.sim,file=paste0("./simulation_inputs/simulate_data_comparison",i,".Rdata"))
}
