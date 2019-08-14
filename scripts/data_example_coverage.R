# ------------------------------------------------------------------------ #
# ------------------- Data analysis: conditional coverage ---------------- #
# ------------------------------------------------------------------------ #

# direc <- "../moretrees/" # path of the moretrees repository
# direc <- "~/Documents/PhD_Papers/Air_pollution/R_code/MORETreeS/moretrees/"
direc <- "/nfs/home/E/ethomas/shared_space/ci3_nsaph/Emma/R_code/moretrees"
setwd(direc)

### Load functions ###
source("scripts/VI_functions.R")
source("scripts/processing_functions.R")
require(igraph)

### Load ICD9 tree ###
load("./simulation_inputs/inputs.Rdata")

######### Algorithm parameters #########

datArgs <- as.integer(as.character(commandArgs(trailingOnly = TRUE))) # Use to call arguments from the command line
# datArgs <- c(0,20,1E5,1E-8,10) # Alternatively, enter arguments directly in R

block <- datArgs[1]+1 # which block of simulations is this
nsims <- datArgs[2] # how many sims to do here
m.max <- datArgs[3] # maximum number of time steps
tol <- datArgs[4] # tolerance for convergence
m.print <- datArgs[5] # how often to print progress for VI algorithm

######################## Load inputs ########################

# Load data
load(file="/nfs/home/E/ethomas/shared_space/ci3_nsaph/Emma/Data/moretrees_data/moretrees_CC_data.Rdata")

# Load original results
load(file="./data_example_results/data_example_full.Rdata")
# Use the collapsed estimates as true betas in simulation
beta_true <- final_ss$moretrees_est
# Use other estimates as initial values for VI algorithm in simulation
mu_gamma_init <- final_ss$VI_params$mu_gamma
sigma2_gamma_init <- final_ss$VI_params$sigma2_gamma
u_s_init <- final_ss$VI_params$u_s
tau_init <- final_ss$hyperparams[2]
rho_init <- final_ss$hyperparams[1]


######################## Variational inference ########################

# For storing results
if(!dir.exists("simulation_results/coverage")) dir.create("simulation_results/coverage")

beta_file <- paste0("simulation_results/coverage/beta_coverage_block",block,".csv")
VI_file <- paste0("simulation_results/coverage/VI_coverage_block",block,".csv")
hyper_file <- paste0("simulation_results/coverage/hyper_coverage_block",block,".csv")

if(file.exists(beta_file)){
  i <- nrow(read.csv(beta_file))
} else {
  file.create(beta_file)
  i <- 0
}
if(!file.exists(VI_file))  file.create(VI_file)
if(!file.exists(hyper_file)) file.create(hyper_file)


# saving output to track number of time steps
pr <- paste0("./simulation_results/coverage/print_coverage_block",block,".txt")
sink(file=pr)

# create list for simulated data
Z.sim <- Z

while(i <= nsims){
  
  # Simulate exposures and case/control status\
  for(v in 1:pL){
    # Compute probability unit 1 is a case
    prob_unit1 <- 1/(1+exp(-beta_true[v]*Z[[v]]))
    # Randomly assign case or control status
    cc <- runif(Y[v]) < prob_unit1
    Z.sim[[v]] <- Z[[v]]*(2*cc-1)
  }
  
  # run VI algorithm
  out_vi <- VI_binary_ss(Z=Z.sim,Y=Y,n=sum(Y),p=p,pL=pL,ancestors=ancestors,
                         leaf.descendants=leaf.descendants,cutoff=0.5,
                         mu_gamma_init=mu_gamma_init,u_s_init=u_s_init,sigma2_gamma_init=sigma2_gamma_init,
                         tau_init=tau_init,rho_init=rho_init,
                         tol=tol,m.max=m.max,m.print=m.print,more=FALSE,update_hyper=T,update_hyper_freq=10)
  
  # iterate sim number
  i <- i + 1
  cat("\n\nCompleted simulation",i,"\n\n")
  
  # write output to csv file
  write.table(rbind(c(i,out_vi$moretrees_est)), file = beta_file, row.names =FALSE, col.names = FALSE,sep = ",", append = TRUE)
  
  write.table(cbind(rep(i,p),out_vi$VI_params), file = VI_file, row.names =FALSE, col.names = FALSE,sep = ",", append = TRUE)
  
  write.table(rbind(c(i,out_vi$hyperparams)), file = hyper_file, row.names =FALSE, col.names = FALSE,sep = ",", append = TRUE)
  
}

sink() # close sink
