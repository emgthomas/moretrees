# --------------------------------------------------------------------- #
# ------------------------- Sensitivity analysis ---------------------- #
# -------- (only one episode of hospitalization per individual) --------#
# --------------------------------------------------------------------- #

# direc <- "../moretrees/" # path of the moretrees repository
# direc <- "/Users/emt380/Documents/PhD_Papers/Air_pollution/R_code/MORETreeS/moretrees/"
direc <- "/nfs/home/E/ethomas/shared_space/ci3_nsaph/Emma/R_code/moretrees"
setwd(direc)

if(!dir.exists("./data_example_results")) dir.create("./data_example_results")

### Load functions ###
source("scripts/VI_functions.R")
source("scripts/processing_functions.R")
require(igraph)
require(doParallel)
require(data.table)

### Load ICD9 tree ###
load("./simulation_inputs/inputs.Rdata")
# Extract list of relevant ICD9 codes
codes <- names(V(tree)[V(tree)$leaf])

######### Algorithm parameters #########

datArgs <- as.integer(as.character(commandArgs(trailingOnly = TRUE))) # Use to call arguments from the command line
# datArgs <- c(3,1E5,1E-16) # Alternatively, enter arguments directly in R

nrestarts <- datArgs[1] # number of random restarts
m.max <- datArgs[2] # maximum number of time steps
tol <- datArgs[3] # tolerance for convergence
m.print <- datArgs[4] # print time step every m.print steps

############### Prepare data ###############

# Load data
load(file="/nfs/home/E/ethomas/shared_space/ci3_nsaph/Emma/Data/moretrees_data/moretrees_CC_data.Rdata")

############ Randomly select one event per enrollee #########
# Get Y, Z, and Q in a format that's easier to work with
Z.unlist <- unlist(Z)
Y.unlist <- unlist(sapply(1:length(Y),function(i,Y) rep(i,Y[i]), Y=Y))
Q.unlist <- unlist(Q)
# # check
# for(i in 1:pL) if(!all.equal(Z[[i]],Z.unlist[Y.unlist==i])) print(paste0("Error in outcome ",i)) else print(i)
# # all good

# Convert to data.table to do the search
dat <- data.table(Q=Q.unlist,Z=Z.unlist,Y=Y.unlist)
# sort by QID
Q.tab <- setkey(dat,Q)
Q.tab$rows <- 1:nrow(Q.tab)
# Set seed
set.seed(7683958)
# Randomly select events to keep (one per enrollee)
Q.keep <- Q.tab[,.I[sample(.N,1)], by=Q]
Q.tab.no.mults <- Q.tab[Q.keep$V1,]
# Sort by outcome (Y)
Q.tab.no.mults <- setkey(Q.tab.no.mults,Y)
# Return to original format
Y.no.mults <- Q.tab.no.mults[,.(.N),by="Y"]$N
# # check
# sum(Y.no.mults <= Y) # yep
Z.no.mults <- Q.tab.no.mults[,list(.(Z)),by="Y"]$V1
# # check
# sum(Y.no.mults == sapply(Z.no.mults,length)) # yep

############### Run analysis on whole dataset ###############

# Adhoc collapsing estimates
adhoc_coeffs <- adhoc_collapsing(Z.no.mults,Y.no.mults,pL,groups)

# Initial values for node coefficients
nodes_init <- initial_node_coeffs(Z.no.mults,Y.no.mults,uncollapsed=adhoc_coeffs[,1],p,pL,leaf.descendants,ancestors)

# For parallelization
registerDoParallel(cores=nrestarts)

# Run spike & slab model using nrestarts random restarts
restarts_ss <- foreach(j = 1:nrestarts) %dopar% {
  
  # saving output to track number of time steps
  pr <- paste0("./data_example_results/data_example_full_sensitivity_print",j,".txt")
  sink(file=pr)
  
  # run VI algorithm
  out_vi <- VI_binary_ss(Z=Z.no.mults,Y=Y.no.mults,n=sum(Y.no.mults),p=p,pL=pL,ancestors=ancestors,
               leaf.descendants=leaf.descendants,cutoff=0.5,mu_gamma_init=nodes_init,
               tol=tol,m.max=m.max,m.print=m.print,more=FALSE,update_hyper=TRUE,update_hyper_freq=10)
  
  sink() # close sink
  
  out_vi
  
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
  Y.g <- sum(Y.no.mults[which.dat])
  if(Y.g == 0){
    beta.groups[which.dat,i] <- NA
  } else {
    beta_ml <- glm(rep(1,Y.g) ~ 0 + unlist(Z.no.mults[which.dat]),family="binomial")
    beta.ml.groups[g,] <- c(beta_ml$coefficients[1],
                          beta_ml$coefficients[1]+qnorm(0.025)*sqrt(vcov(beta_ml)),
                          beta_ml$coefficients[1]+qnorm(0.975)*sqrt(vcov(beta_ml)))
  }
}

############### Save results ###############

save(beta_est,groups,beta.ml.groups,Y.no.mults,final_ss,adhoc_coeffs,ELBOS,file = "data_example_results/data_example_full_sensitivity.Rdata")
