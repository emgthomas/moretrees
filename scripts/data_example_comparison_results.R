# ------------------------------------------------------------------------ #
# ------------------- Results for VI/MCMC comparison --------------------- #
# ------------------------------------------------------------------------ #

# direc <- "/nfs/home/E/ethomas/shared_space/ci3_nsaph/Emma/R_code/moretrees"
direc <- "/Users/emt380/Documents/PhD_Papers/Air_pollution/R_code/MORETreeS/moretrees"
# direc <- "../moretrees/" # path of the moretrees repository
setwd(direc)

#### Create directory for saving results ###
if(!dir.exists("./figures_and_tables")) dir.create("./figures_and_tables")

### Load functions
require(bayesplot)
require(igraph)
require(Matrix)
require(BoomSpikeSlab)
source("scripts/processing_functions.R")

### Load ICD9 tree
load("./simulation_inputs/inputs.Rdata")

# Get ancestor matrix A
A <- t(as_adj(tree,sparse = T))
A <- expm(A)
A[A>0] <- 1 
A.conv <- A[(p-pL+1):p,] # matrix for conversion: beta <- A.conv %*% gamma

############################### Convergence diagnostics ###############################

nsims <- 10
nchains <- 3
beta_mcmc <- list()
gamma_mcmc <- list()
for(sim in 1:nsims){
  # beta <- list()
  # gamma <- list()
  beta <- Matrix(nrow=0,ncol=pL+1)
  gamma <- Matrix(nrow=0,ncol=p+1)
  for(chain in 1:nchains){ 
    # if(!(sim==4 & chain==1)){ # remove this later
    # load mcmc results
    res_mcmc <- paste0("./data_example_results/comparison_mcmc_samples",sim,"_chain",chain,".csv")
    iters <- Matrix(as.matrix(read.csv(file=res_mcmc,header = F)))
    # if(sim ==4) chain <- chain-1
    # gamma[[chain]] <- mcmc(iters)
    gamma <- rbind(gamma,cbind(iters,rep(chain,nrow(iters))))
    beta_chain <- iters %*% t(A.conv)
    beta <- rbind(beta,cbind(beta_chain,rep(chain,nrow(beta_chain))))
    # beta[[chain]] <- mcmc(beta_chain)
    # }
  }
  beta.df <- as.data.frame(as.matrix(beta))
  names(beta.df)[pL+1] <- "Chain"
  beta_mcmc[[sim]] <- beta.df
  gamma.df <- as.data.frame(as.matrix(gamma))
  names(gamma.df)[1:p] <- colnames(A)
  names(gamma.df)[p+1] <- "Chain"
  gamma_mcmc[[sim]] <- gamma.df
  # beta_mcmc[[sim]] <- mcmc.list(beta)
  # gamma_mcmc[[sim]] <- mcmc.list(gamma)
  print(sim)
}

sim <- 9
u <- c("44100","44101","44102","44103")
mcmc_trace(beta_mcmc[[sim]],par=u)
mcmc_pairs(beta_mcmc[[sim]],par=u)
# u <- c(200,201)
# mcmc_trace(beta_mcmc[[sim]],par=names(beta_mcmc[[sim]])[u])
# mcmc_pairs(beta_mcmc[[sim]],par=names(beta_mcmc[[sim]])[u])

v <- c("44100","44101","44102","44103")
# mcmc_trace(gamma_mcmc[[sim]],par=names(gamma_mcmc[[sim]])[v]) 
# mcmc_pairs(gamma_mcmc[[sim]],par=names(gamma_mcmc[[sim]])[v])
mcmc_trace(gamma_mcmc[[sim]],par=v) 
mcmc_pairs(gamma_mcmc[[sim]],par=v)


# require(coda)
# gelman.diag(beta_mcmc[[1]])

############################### Mean and variance comparison ###############################

var_VI <- matrix(nrow=nsims,ncol=pL)
var_MCMC <- matrix(nrow=nsims,ncol=pL)
for(sim in 1:nsims){
  # load VI results
  res_VI <- paste0("./data_example_results/comparison_vi_results",sim,".Rdata")
  load(res_VI)
  out_vi$VI_params$sim <- 1
  beta_sd <- sapply(1:p,FUN=indiv.beta.sd.calc,VIsims=out_vi$VI_params,ancestors=ancestors,pL=pL,p=p)
  var_VI[sim,] <- 
}
