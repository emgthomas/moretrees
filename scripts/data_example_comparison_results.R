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
require(glue)
require(RColorBrewer)
require(patchwork)
require(icd)
# # To install the patchwork package:
# library(devtools)
# install_github("thomasp85/patchwork")
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
  names(beta.df)[1:pL] <- colnames(A)[(p-pL+1):p]
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

#### Trace plots
for(sim in 2:nsims){
  plots.list <- list()
  codes <- names(beta_mcmc[[sim]])[1:pL]
  for(u in 1:pL){
    outcome <- codes[u]
    plots.list[[u]] <- mcmc_trace(beta_mcmc[[sim]],
                                  par=outcome,
                                  size=0.2) + 
      theme_minimal() +
      theme(legend.position="none",
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.grid = element_blank(),
            plot.title=element_text(hjust=0.5,size=10)) +
      ggtitle(as.character(icd_short_to_decimal(outcome)))
  }
  # reordering
  plots.list <- plots.list[order(codes)]
  # pdf(file=paste0("./figures_and_tables/figureA6_",sim,".pdf"),width=16,height=27)
  png(file=paste0("./figures_and_tables/figureA6_",sim,".png"),width=16,height=27,units="in",
      res=500)
  print(wrap_plots(plots.list,ncol=16))
  dev.off()
}

############################### Mean and variance comparison ###############################

### compute beta estimates ###
beta_est_mcmc <- data.frame(matrix(nrow=pL,ncol=nsims))
names(beta_est_mcmc) <- names(A)[(p-pL+1):p]
group_est_mcmc <- data.frame(matrix(nrow=pL,ncol=nsims))
names(group_est_mcmc) <- names(A)[(p-pL+1):p]
beta_est_vi <- data.frame(beta=matrix(nrow=pL,ncol=nsims))
names(beta_est_vi) <- names(A)[(p-pL+1):p]
group_est_vi <- data.frame(beta=matrix(nrow=pL,ncol=nsims))
names(group_est_vi) <- names(A)[(p-pL+1):p]
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
  group_est_mcmc[,sim] <- as.integer(factor(beta_est,
                                          levels=sort(unique(beta_est))))
  ## VI MOReTreeS estimate
  res_VI <- paste0("./data_example_results/comparison_vi_results",sim,".Rdata")
  load(res_VI)
  beta_est_vi[,sim] <- out_vi$moretrees_est
  group_est_vi[,sim] <- as.integer(factor(out_vi$moretrees_est,
                                          levels=sort(unique(out_vi$moretrees_est))))
  p_s_vi <- 1/(1+exp(-out_vi$VI_params$u_s))
}

########## Compare VI and MCMC estimates of beta ##############

### Load data example results
load("./data_example_results/data_example_full.Rdata")

### Sample size
n <- 1E6 # Size of simulated dataset
n_data <- sum(Y) # Size of original dataset
m <- sqrt(n_data/n) # Scaling factor for coefficients

### Beta estimates used for simultaions
beta_indiv_est <- indiv.beta.ci.calc(final_ss$VI_params,ancestors,pL,p)
# Scale up the beta estimates to account for smaller sample size
beta_sim <- beta_indiv_est$beta_indiv*m

####### Prepare data for plotting ########
vi.mcmc.smc <- matrix(nrow=0,ncol=3)
mcmc.vi.smc <- matrix(nrow=0,ncol=3)

vi.dat <- data.frame(est=beta_est_vi,group=group_est_vi,method="VI")
mcmc.dat <- data.frame(est=beta_est_mcmc,group=group_est_mcmc,method="MCMC")
dat <- rbind(vi.dat,mcmc.dat)
dat <- data.frame(est_vi=beta_est_vi,est_mcmc=beta_est_mcmc,
                  group_vi=group_est_vi,group_mcmc=group_est_mcmc)
dat.df <- reshape(dat,direction="long",
                  varying=list(1:10,11:20,21:30,31:40))
names(dat.df)[1:6] <- c("sim","est_vi","est_mcmc","group_vi","group_mcmc","node")
dat.df$est_true <- beta_sim
change.group <- apply(X=cbind(as.character(dat.df$sim-1),
                              as.character(dat.df$group_vi),
                              as.character(dat.df$group_mcmc)),
                      MARGIN=1,FUN=glue_collapse,sep="")
dat.df$change.group <- change.group
n.outcomes <- tapply(change.group,as.factor(change.group),length)
# n.cases <- tapply(dat.df$Y,as.factor(change.group),sum)
n.df <- data.frame(change.group=names(n.outcomes),
                   n.outcomes=n.outcomes)
                   #, n.cases=n.cases)
dat.df <- merge(dat.df,n.df,by="change.group",all.x=T,all.y=F)
dat.df$est_vi <- exp(dat.df$est_vi)
dat.df$est_mcmc <- exp(dat.df$est_mcmc)
dat.df$est_vi_lab <- sprintf("%.4f",dat.df$est_vi)
dat.df$est_mcmc_lab <- sprintf("%.4f",dat.df$est_mcmc)
dat.df$est_true_lab <- sprintf("%.4f",dat.df$est_true)

##### Plotting ORs estimated by MCMC vs those estimated by VI
top.pts <- 20
left.pts <- 20
plot.list <- list()
for(sim in 1:10){
  # groups plot
  dat.sim <- dat.df[dat.df$sim==sim,]
  dat.sim$est_vi_lab <- factor(dat.sim$est_vi_lab,
                           levels=sort(unique(dat.sim$est_vi_lab),decreasing = F))
  dat.sim$est_mcmc_lab <- factor(dat.sim$est_mcmc_lab,
                           levels=sort(unique(dat.sim$est_mcmc_lab),decreasing = F))
  plot.groups <- ggplot(dat.sim,
                        aes(x=est_mcmc_lab,y=est_vi_lab,label=n.outcomes)) + 
    geom_point(color="white",size=6) +
    geom_label(size=4,label.size=0,label.padding=unit(0,"lines")) +
    theme_bw() +
    theme(legend.position="none",
          plot.margin = margin(t=-15,r=0,b=0,l=left.pts, unit = "pt"),
          # plot.margin = margin(t=0,r=0,b=10,l=10, unit = "pt"),
          plot.title = element_blank(),
          axis.text=element_text(size=7)) +
    xlab("MCMC Estimates") +
    ylab("VI Estimates")
  
  # simple matching coefficient plot
  ### VI
  smc.vi.df <- data.frame(OR=sort(unique(dat.sim$est_vi_lab)),
                          smc=smc(dat.sim$group_vi,dat.sim$group_mcmc))
  smc.vi.df$label <- sprintf("%.3f",smc.vi.df$smc)
  smc.vi.df$label[smc.vi.df$label=="1"] <- rep("1.000",nrow(smc.vi.df))
  plot.vi.smc <- ggplot(smc.vi.df,aes(y=OR,x=1)) + 
    geom_tile(aes(fill=smc),colour="black",size=0.1) +
    geom_label(aes(x=1,y=OR,label=label),size=2,alpha=0.6,label.size=0,label.padding=unit(0.1,"lines")) +
    theme_bw() +
    theme(legend.position="none",
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.margin = margin(t=-15,r=0,b=0,l=0, unit = "pt"),
          # plot.margin = margin(t=-3,r=0,b=5,l=-3, unit = "pt"),
          plot.title = element_blank(),
          panel.border=element_blank()) +
    xlab("") + 
    scale_fill_gradientn(colours=brewer.pal(7,name="YlGnBu"),
                         limits=c(0,1)) +
    scale_x_continuous(breaks=0.5,labels="")
  
  ### mcmc
  smc.mcmc.df <- data.frame(OR=sort(unique(dat.sim$est_mcmc_lab)),
                          smc=smc(dat.sim$group_mcmc,dat.sim$group_vi))
  smc.mcmc.df$label <- sprintf("%.3f",smc.mcmc.df$smc)
  smc.mcmc.df$label[smc.mcmc.df$label=="1"] <- rep("1.000",nrow(smc.mcmc.df))
  plot.mcmc.smc <- ggplot(smc.mcmc.df,aes(y=1,x=OR)) + 
    geom_tile(aes(fill=smc),colour="black",size=0.1) +
    geom_label(aes(y=1,x=OR,label=label),
               size=2,alpha=0.6,label.size=0,label.padding=unit(0.1,"lines")) +
    theme_bw() +
    theme(legend.position="none",
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          # plot.margin = margin(t=0,r=-4,b=-16,l=52, unit = "pt"),
          plot.margin = margin(t=top.pts,r=0,b=0,l=left.pts, unit = "pt"),
          panel.border=element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    xlab("") + 
    ylab("") +
    scale_fill_gradientn(colours=brewer.pal(7,
                                            name="YlGnBu"),
                         limits=c(0,1)) +
    ggtitle(paste0("Simulation ",sim))
  
  # filler plot for corner
  dat.fill <- data.frame(x=1,y=1,label="kappa")
  plot.filler <- ggplot(dat.fill,aes(x=x,y=y,label=label)) +
    theme_bw() +
    geom_text(parse = TRUE,size=6)+
    theme(legend.position="none",
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.margin = margin(t=top.pts,r=0,b=0,l=0, unit = "pt"),
          panel.border=element_blank())
  
  # Save plots as pdf
  plot.list[[sim]] <- plot.mcmc.smc + 
    plot.filler + 
    plot.groups + 
    plot.vi.smc +
    plot_layout(ncol=2,widths=c(7,1),heights=c(1,7))

  print(sim)
}

pdf(file = paste0("./figures_and_tables/figureA5.pdf"),width=20,height=8)
wrap_plots(plot.list,ncol=5,widths=rep(1,5),heights=rep(1,2))
dev.off()

### compute variance estimates ###
nonzero_var <- function(x){
  if(sum(x!=0)<=1) return(0)
  return(var(x[x!=0]))
} 
var_VI <- matrix(nrow=pL,ncol=nsims)
var_MCMC <- matrix(nrow=pL,ncol=nsims)
var_nonzero_mat <- matrix(nrow=0,ncol=4)
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
  var_nonzero_mat <- rbind(var_nonzero_mat,
                         cbind(var_VI_nonzero,
                               var_MCMC_nonzero,
                               which_nonzero,
                               rep(sim,length(var_VI_nonzero))))
  row.names(var_nonzero_mat) <- NULL
}

### Plot variance for all components
var_dat <- data.frame(MCMC=var_MCMC,VI=var_VI,node=1:pL)
var_dat <- reshape(var_dat,varying=list(1:10,11:20),direction="long")
var_dat$id <- NULL
var_dat$SmallerVariance <- "MCMC"
var_dat$SmallerVariance[var_dat$VI <= var_dat$MCMC] <- "VI"
names(var_dat) <- c("node","sim","MCMC","VI","SmallerVariance")

axis.min <- min(as.matrix(var_dat[,c("VI","MCMC")]))
axis.max <- max(as.matrix(var_dat[,c("VI","MCMC")]))
line.x <- c(axis.min,axis.max)
var_plot <- ggplot(var_dat,aes(x=VI,y=MCMC,col=SmallerVariance)) + 
  geom_abline(intercept=0,slope=1) +
  geom_point(alpha=0.8,shape=1) +
  facet_wrap(.~sim,ncol=5) +
  scale_x_log10(lim=c(axis.min,axis.max),breaks=c(1E-6,1)) +
  scale_y_log10(lim=c(axis.min,axis.max),breaks=c(1E-6,1)) +
  theme_minimal()

pdf(file="./figures_and_tables/figureA7.pdf",width=10,height=4)
var_plot
dev.off()

### Plot variance for non-zero components
var_nonzero_dat <- as.data.frame(var_nonzero_mat)
names(var_nonzero_dat) <- c("VI","MCMC","nodes","sim")
var_nonzero_dat$SmallerVariance <- "MCMC"
var_nonzero_dat$SmallerVariance[var_nonzero_dat$VI <= var_nonzero_dat$MCMC] <- "VI"

axis.min <- min(as.matrix(var_nonzero_dat[,c("VI","MCMC")]))
axis.max <- max(as.matrix(var_nonzero_dat[,c("VI","MCMC")]))
line.x <- c(axis.min,axis.max)
var_nonzero_plot <- ggplot(var_nonzero_dat,aes(x=VI,y=MCMC,col=SmallerVariance)) + 
  geom_abline(intercept=0,slope=1) +
  geom_point(alpha=1,shape=1,size=3) +
  facet_wrap(.~sim,ncol=5) +
  scale_x_log10(lim=c(axis.min,axis.max),breaks=c(1E-7,1E-5)) + 
  scale_y_log10(lim=c(axis.min,axis.max),breaks=c(1E-7,1E-5)) +
  theme_minimal()

pdf(file="./figures_and_tables/figureA8.pdf",width=10,height=4)
var_nonzero_plot
dev.off()

### Average variance ratio VI:MCMC across sims for nonzero components
var_nonzero_dat$ratio <- var_nonzero_dat$VI/var_nonzero_dat$MCMC
ratio_sims <- tapply(var_nonzero_dat$ratio,
                     var_nonzero_dat$sim,
                     mean)
cat(paste0("\n\nOn average, VI estimate of variance for nonzero gammas was ",
           format(100*(1-mean(ratio_sims)),digits=3),
           "% lower than MCMC estimate"))


