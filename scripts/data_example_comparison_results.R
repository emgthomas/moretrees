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
require(reshape2)
require(ggplot2)
source("scripts/processing_functions.R")

### Load ICD9 tree
load("./simulation_inputs/inputs.Rdata")

# Get ancestor matrix A
A <- t(as_adj(tree,sparse = T))
A <- expm(A)
A[A>0] <- 1 
A.conv <- A[(p-pL+1):p,] # matrix for conversion: beta <- A.conv %*% gamma

############################### Loading MCMC chains ###############################

nsims <- 10
nchains <- 3
beta_mcmc <- list()
gamma_mcmc <- list()
p_nonzero_mcmc <- list()
for(sim in 1:nsims){
  beta <- Matrix(nrow=0,ncol=pL+1)
  gamma <- Matrix(nrow=0,ncol=p+1)
  p_nonzero <- Matrix(nrow=0,ncol=p+1)
  for(chain in 1:nchains){ 
    res_mcmc <- paste0("./data_example_results/comparison_mcmc_samples",sim,"_chain",chain,".csv")
    iters <- Matrix(as.matrix(read.csv(file=res_mcmc,header = F)))
    gamma <- rbind(gamma,cbind(iters,rep(chain,nrow(iters))))
    beta_chain <- iters %*% t(A.conv)
    beta <- rbind(beta,cbind(beta_chain,rep(chain,nrow(beta_chain))))
    p_nonzero <- rbind(p_nonzero,cbind(iters!=0,rep(chain,nrow(iters))))
  }
  beta.df <- as.data.frame(as.matrix(beta))
  names(beta.df)[1:pL] <- names(Z.conv)
  names(beta.df)[pL+1] <- "Chain"
  beta_mcmc[[sim]] <- beta.df
  gamma.df <- as.data.frame(as.matrix(gamma))
  names(gamma.df)[1:p] <- colnames(A)
  names(gamma.df)[p+1] <- "Chain"
  gamma_mcmc[[sim]] <- gamma.df
  p_nonzero.df <- as.data.frame(as.matrix(p_nonzero))
  names(p_nonzero.df)[1:p] <- colnames(A)
  names(p_nonzero.df)[p+1] <- "Chain"
  p_nonzero_mcmc[[sim]] <- p_nonzero.df
  print(sim)
}

############################### Convergence diagnostics ###############################

sim <- 9
u <- c("44100","44101","44102","44103")
mcmc_trace(beta_mcmc[[sim]],par=u)
mcmc_pairs(beta_mcmc[[sim]],par=u)
# u <- c(200,201)
# mcmc_trace(beta_mcmc[[sim]],par=names(beta_mcmc[[sim]])[u])
# mcmc_pairs(beta_mcmc[[sim]],par=names(beta_mcmc[[sim]])[u])

v <- c("44100","44101","44102","44103")
mcmc_trace(gamma_mcmc[[sim]],par=v) 
mcmc_pairs(gamma_mcmc[[sim]],par=v)
v <- 1
# mcmc_trace(gamma_mcmc[[sim]],par=names(gamma_mcmc[[sim]])[v]) 
# mcmc_pairs(gamma_mcmc[[sim]],par=names(gamma_mcmc[[sim]])[v])

############################### Mean and variance comparison ###############################

### compute beta estimates ###
beta_est_mcmc <- data.frame(matrix(nrow=pL,ncol=nsims))
names(beta_est_mcmc) <- names(A.conv)
beta_est_vi <- data.frame(beta=matrix(nrow=pL,ncol=nsims))
names(beta_est_vi) <- names(A.conv)
cutoff <- 0.5
p_s <- matrix(nrow=p,ncol=nsims)
for(sim in 1:nsims){
  ## MCMC MOReTreeS estimate
  p_s[,sim] <- colMeans(p_nonzero_mcmc[[sim]][,1:p])
  nonzero <- p_s[,sim] >= cutoff
  sgamma <- apply(gamma_mcmc[[sim]][,1:p],2,function(x) mean(x[x!=0]))
  sgamma[!nonzero] <- 0
  beta_est <- as.numeric(sgamma %*% t(A.conv))
  beta_est_mcmc[,sim] <- beta_est
  
  ## VI MOReTreeS estimate
  res_VI <- paste0("./data_example_results/comparison_vi_results",sim,".Rdata")
  load(res_VI)
  beta_est_vi[,sim] <- out_vi$moretrees_est
  p_s_vi <- 1/(1+exp(-out_vi$VI_params$u_s))
}

### compute variance estimates ###
nonzero_var <- function(x){
  if(sum(x!=0)<=1) return(0)
  return(var(x[x!=0]))
} 
var_VI <- matrix(nrow=pL,ncol=nsims)
var_MCMC <- matrix(nrow=pL,ncol=nsims)
var_plots_dat <- matrix(nrow=0,ncol=4)
for(sim in 1:nsims){
  # load VI results
  res_VI <- paste0("./data_example_results/comparison_vi_results",sim,".Rdata")
  load(res_VI)
  VI_params <- out_vi$VI_params
  # Variance of betas
  var_VI[,sim] <- beta.var.calc(mu_gamma=VI_params$mu_gamma,
                                sigma2_gamma=VI_params$sigma2_gamma,
                                u_s=VI_params$u_s,
                                ancestors=ancestors,
                                pL=pL,p=p)
  var_MCMC[,sim] <- apply(beta_mcmc[[sim]][,1:pL],2,FUN=var)
  # Variance of non-zero gammas
  which_nonzero <- which(p_s[,sim]>=cutoff)
  var_VI_nonzero <- VI_params$sigma2_gamma[which_nonzero]
  var_MCMC_nonzero <- apply(gamma_mcmc[[sim]][,1:p],2,FUN=nonzero_var)[which_nonzero]
  var_plots_dat <- rbind(var_plots_dat,
                         cbind(var_VI_nonzero,
                               var_MCMC_nonzero,
                               which_nonzero,
                               rep(sim,length(var_VI_nonzero))))
  row.names(var_plots_dat) <- NULL
}

### Plot variance for non-zero components
var_plots_dat <- as.data.frame(var_plots_dat)
names(var_plots_dat) <- c("VI","MCMC","nodes","sim")
var_plots_dat$VI_smaller <- var_plots_dat$VI <= var_plots_dat$MCMC

axis.min <- 0
axis.max <- max(as.matrix(var_plots_dat[,c("VI","MCMC")]))
line.x <- c(axis.min,axis.max)
var_plot <- ggplot(var_plots_dat,aes(x=VI,y=MCMC,col=VI_smaller)) + 
  geom_abline(intercept=0,slope=1) +
  geom_point() +
  facet_wrap(.~sim,ncol=5) +
  scale_x_continuous(lim=c(axis.min,axis.max)) +
  scale_y_continuous(lim=c(axis.min,axis.max))


