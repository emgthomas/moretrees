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
# datArgs <- c(0,2,3,1E-8,1) # Alternatively, enter arguments directly in R

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

betaML_file <- paste0("simulation_results/coverage/beta_ML_coverage_block",block,".csv")
VI_file <- paste0("simulation_results/coverage/VI_coverage_block",block,".csv")
hyper_file <- paste0("simulation_results/coverage/hyper_coverage_block",block,".csv")

if(file.exists(betaML_file) & file.exists(VI_file) & file.exists(hyper_file)){
  i <- nrow(read.csv(hyper_file,header=F))
# } else {
  file.create(betaML_file)
  file.create(VI_file)
  file.create(hyper_file)
  i <- 0
}


# saving output to track number of time steps
pr <- paste0("./simulation_results/coverage/print_coverage_block",block,".txt")
sink(file=pr)

# create list for simulated data
Z.sim <- Z

while(i <= nsims){
  
  # Simulate exposures and case/control status
  for(v in 1:pL){
    # Compute probability unit 1 is a case
    prob_unit1 <- 1/(1+exp(-beta_true[v]*Z[[v]]))
    # Randomly assign case or control status
    cc <- runif(Y[v]) < prob_unit1
    Z.sim[[v]] <- Z[[v]]*(2*cc-1)
  }
  
  # add a little perturbation to initial values
  f <- 0.1
  mu_gamma_init2 <- mu_gamma_init + sapply(abs(mu_gamma_init)*f,rnorm,mean=0,n=1)
  u_s_init2 <- u_s_init + sapply(abs(u_s_init)*f,rnorm,mean=0,n=1)
  sigma2_gamma_init2 <- sigma2_gamma_init + sapply(abs(sigma2_gamma_init)*f,rnorm,mean=0,n=1)
  tau_init2 <- exp(log(tau_init) + rnorm(1,sd=0.1*abs(log(tau_init))))
  rho_init2 <- rho_init + runif(1,rho_init/2,rho_init*3/2)
  
  # run VI algorithm
  out_vi <- VI_binary_ss(Z=Z.sim,Y=Y,n=sum(Y),p=p,pL=pL,ancestors=ancestors,
                         leaf.descendants=leaf.descendants,cutoff=0.5,
                         mu_gamma_init=mu_gamma_init2,u_s_init=u_s_init2,sigma2_gamma_init=sigma2_gamma_init2,
                         tau_init=tau_init2,rho_init=rho_init2,
                         tol=tol,m.max=m.max,m.print=m.print,more=FALSE,update_hyper=T,update_hyper_freq=10)
  
  # fit MLE to discovered groups
  beta_est <- out_vi$moretrees_est
  groups <- as.numeric(as.factor(beta_est))
  beta.ml.groups <- numeric(pL)
  beta.cil.groups <- numeric(pL)
  beta.ciu.groups <- numeric(pL)
  for(g in 1:max(groups)){
    which.dat <- groups==g
    Y.g <- sum(Y[which.dat])
    if(Y.g == 0){
      beta.groups[which.dat,i] <- NA
    } else {
      beta_ml <- glm(rep(1,Y.g) ~ 0 + unlist(Z[which.dat]),family="binomial")
      beta.ml.groups[which.dat] <- beta_ml$coefficients[1]
      beta.cil.groups[which.dat] <- beta_ml$coefficients[1]+qnorm(0.025)*sqrt(vcov(beta_ml))
      beta.ciu.groups[which.dat] <- beta_ml$coefficients[1]+qnorm(0.975)*sqrt(vcov(beta_ml))
    }
  }
  
  # iterate sim number
  i <- i + 1
  cat("\n\nCompleted simulation",i,"\n\n")
  
  # write output to csv file
  write.table(cbind(rep(i,pL),beta.ml.groups,beta.cil.groups,beta.ciu.groups), file = betaML_file, row.names =FALSE, col.names = FALSE,sep = ",", append = TRUE)
  
  write.table(cbind(rep(i,p),out_vi$VI_params), file = VI_file, row.names =FALSE, col.names = FALSE,sep = ",", append = TRUE)
  
  write.table(rbind(c(i,out_vi$hyperparams)), file = hyper_file, row.names =FALSE, col.names = FALSE,sep = ",", append = TRUE)
  
}

sink() # close sink
