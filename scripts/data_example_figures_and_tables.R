# --------------------------------------------------------------------------------------------- #
# ----------------------- Producing figures and tables for data example ----------------------- #
# --------------------------------------------------------------------------------------------- #

# direc <- "../moretrees/" # path of the moretrees repository
# direc <- "/nfs/home/E/ethomas/shared_space/ci3_nsaph/Emma/R_code/moretrees/"
direc <- "/Users/emt380/Documents/PhD_Papers/Air_pollution/R_code/MORETreeS/moretrees/"
setwd(direc)

#### Create directory for saving results ###
if(!dir.exists("./figures_and_tables")) dir.create("./figures_and_tables")

# Load functions
source("./scripts/VI_functions.R")
source("./scripts/processing_functions.R")
require(Matrix)
require(igraph)
require(ggplot2)
require(igraph)
require(mclust)
require(patchwork)
require(data.table)
require(reshape2)
require(icd)
require(glue)
require(RColorBrewer)
require(collapsibleTree)
require(circlize)
require(plotly)
require(xtable)

############################### Comparison of tree and MLEs ###############################
# Load tree
load("simulation_inputs/inputs.Rdata")
# Load full data results
load("data_example_results/data_example_full.Rdata")

# Compute prior correlation based on tree
outcome_pairs <- data.frame(t(combn((p-pL+1):p,m=2)))
names(outcome_pairs) <- c("u","v")
outcome_pairs$prior_corr <- 
# outcome_pairs$prior_corr <- mapply(FUN=prior_corr_betas,u=outcome_pairs$u,outcome_pairs$v,MoreArgs=list(ancestors=ancestors))
# outcome_pairs$prior_corr_rounded <- round(outcome_pairs$prior_corr,digits=1)
# outcome_pairs$prior_corr_rounded <- factor(outcome_pairs$prior_corr_rounded,
#                                               levels=seq(0.2,0.8,0.1))

# Compute difference in individual MLE estimates
outcome_pairs$mle_u <- adhoc_coeffs$uncollapsed[outcome_pairs$u-(p-pL)]
outcome_pairs$mle_v <- adhoc_coeffs$uncollapsed[outcome_pairs$v-(p-pL)]
outcome_pairs$mle_diff <- abs(outcome_pairs$mle_u-outcome_pairs$mle_v)

# Sample size
outcome_pairs$n_u <- Y[outcome_pairs$u]
outcome_pairs$n_v <- Y[outcome_pairs$v]

# Box plots
prior_plot <- ggplot(outcome_pairs,
                     aes(group=as.factor(prior_corr),
                         x=prior_corr,y=mle_diff)) + 
  geom_smooth(aes(x=prior_corr,y=mle_diff,group=NA),method = "lm",
              colour="blue",lty=2,lwd=0.5,se=TRUE) +
  geom_boxplot(outlier.size=1,outlier.alpha=0.3,alpha=0.5) +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  xlab(bquote("Corr(" ~ beta[u] ~","~ beta[v] ~ ")")) + 
  ylab(bquote("|" ~ beta[v]^{MLE} - beta[u]^{MLE} ~ "|"))

pdf(file="./figures_and_tables/figure4a.pdf",width=5,height=3)
prior_plot
dev.off()

# Consider subset of outcomes with minimum sample size
n_thresh <- 10000
outcome_pairs_subset <- subset(outcome_pairs,n_u >= n_thresh & n_v >= n_thresh)
sum(Y>= n_thresh)

# Box plots
prior_plot2 <- ggplot(outcome_pairs_subset,
                      aes(group=as.factor(prior_corr),
                          x=prior_corr,y=mle_diff))+ 
  geom_smooth(aes(x=prior_corr,y=mle_diff,group=NA),method = "lm",
              colour="blue",lty=2,lwd=0.5,se=TRUE) +
  geom_boxplot(outlier.size=1,outlier.alpha=0.5,alpha=0.5) +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  xlab(bquote("Corr(" ~ beta[u] ~","~ beta[v] ~ ")")) + 
  ylab(bquote("|" ~ beta[v]^{MLE} - beta[u]^{MLE} ~ "|"))

pdf(file="./figures_and_tables/figure4b.pdf",width=5,height=3)
prior_plot2
dev.off()


############################### Cross-validation results ###############################

# Load tree
load("simulation_inputs/inputs.Rdata")

### Select value of tuning parameter ###

nfolds <- 10
load(paste0("data_example_results/data_example_cv_fold1.Rdata"))
ll.tp.folds <- as.data.frame(matrix(nrow=length(tp),ncol=nfolds))
for(i in 1:nfolds){
  load(paste0("data_example_results/data_example_cv_fold",i,".Rdata"))
  ll.tp.folds[,i] <- ll.tp
  #plot(tp,ll.tp,type="l",main = i)
}

# Get average LL over folds for each value of tuning parameter
ll.tp <- rowMeans(ll.tp.folds)

# Best value of tuning parameter
ll.max <- max(ll.tp)
tp.max <- tp[which(ll.tp==ll.max)]
tp.max.mid <- tp.max[round(mean(1:length(tp.max)))]
ll.max.folds <- ll.tp.folds[which(tp==tp.max.mid),]
ll.5.folds <- ll.tp.folds[which(tp==0.5),]

### How many groups result from each value of tp, in full data? ###

# Load full data results
load("data_example_results/data_example_full.Rdata")

# Get ancestor matrix A
A <- t(as_adj(tree,sparse = T))
A <- expm(A)
A[A>0] <- 1
A <- A[(p-pL+1):p,]

# Extract relevant VI params
mu_gamma <- final_ss$VI_params$mu_gamma
p_vi <- exp(loglogit(final_ss$VI_params$u_s))

# Get number of groups for each value of tp
ngroups <- numeric(length(tp))
for(i in 1:length(tp)){
  sgamma <- (mu_gamma * (p_vi >= tp[i]))
  beta.tp <- A %*% sgamma
  ngroups[i] <- length(unique(as.numeric(beta.tp)))
}


# Best beta
sgamma <- (mu_gamma * (p_vi >= min(tp.max)))
beta.tp.min <- A %*% sgamma
sgamma <- (mu_gamma * (p_vi >= max(tp.max)))
beta.tp.max <- A %*% sgamma

# Is it same as beta_coll?
all.equal(as.numeric(beta.tp.min),final_ss$moretrees_est)
all.equal(as.numeric(beta.tp.max),final_ss$moretrees_est)

# Best number of groups
ngroups_max <- unique(ngroups[which(ll.tp==ll.max)])

# Some plotting parameters
grid_lwd <- 0.2
fontsize <- 19
tp.marks <- c(tp[1],0.122,min(tp.max),max(tp.max),0.999)
ll.marks <- c(ll.tp[1],ll.max,ll.tp[length(ll.tp)])
dat <- data.frame(tp=tp,ngroups=ngroups,ll.tp=ll.tp)

### Plot tuning parameter versus CV ll
plot_tp <- ggplot(dat,aes(x=tp,y=ll.tp)) +
  # Plot line
  geom_line() +
  # Aesthetics
  theme_bw(base_size=fontsize) +
  # Grid lines
  theme(panel.grid.major=element_line(size=0.5,linetype=2),
        panel.grid.minor=element_blank()) + 
  # Axes
  scale_x_continuous(limits=c(0,1),breaks=tp.marks) +
  scale_y_continuous(limits=c(min(ll.tp),max(ll.tp)),
                     breaks=ll.marks) +
  xlab(expression(delta)) +
  ylab("Cross-validated mean log likelihood")

### Plot number of groups versus CV ll
ngroups.marks <- c(ngroups[1],7,ngroups_max,26,ngroups[length(ngroups)])
plot_ngroups <- ggplot(dat,aes(x=ngroups,y=ll.tp)) +
  # Plot line
  geom_line() +
  # Aesthetics
  theme_bw(base_size=fontsize) +
  # Grid lines
  theme(panel.grid.major=element_line(size=0.5,linetype=2),
        panel.grid.minor=element_blank()) + 
  # Axes
  scale_x_log10(breaks=ngroups.marks) +
  scale_y_continuous(limits=c(min(ll.tp),max(ll.tp)),
                     breaks=ll.marks) +
  xlab(expression("G("*delta*")")) +
  ylab("Cross-validated mean log likelihood")

############### Figure R.1 (reviewer response) ################
pdf("figures_and_tables/figure_cv_delta_a.pdf",width=8,height=5)
plot_tp
dev.off()

pdf("figures_and_tables/figure_cv_delta_b.pdf",width=8,height=5)
plot_ngroups
dev.off()

### Load other CV results ###
nmods <- 8
cv.res <- as.data.frame(matrix(nrow=nfolds,ncol=nmods+1))
names(cv.res) <- c("fold",
                   "ssMOReTreeS\n Collapsed",
                   #"ssMOReTreeS\n Collapsed\n Max Delta",
                   "ssMOReTreeS\n Individual",
                   "Uncollapsed",
                   "sim_groups",
                   "Adhoc\n Grouping 1",
                   "Adhoc\n Grouping 2",
                   "Adhoc\n Grouping 3",
                   "Fully\n Collapsed")
cv.moretrees <- numeric(nfolds)
for(i in 1:nfolds){
  load(paste0("data_example_results/data_example_cv_fold",i,".Rdata"))
  # cv.res[i,] <- c(ll.cv[1,1:2],ll.max.folds[1,i],as.numeric(ll.cv[1,3:length(ll.cv)]))
  cv.res[i,] <- ll.cv
}
cv.res$sim_groups <- NULL
cv.df <- melt(cv.res,id.vars="fold",measure.vars=2:(nmods),variable.name="mod",value.name="ll")

# Cross validated predictive accuracy
sort(colMeans(cv.res[,2:ncol(cv.res)]),decreasing = T)

# Winner excluding MOReTreeS individual
apply(cv.res[,-c(1,3)],1,which.max)
# How often does MOReTreeS individual beat MOReTrees collapsed?
sum(cv.res[,2] <= cv.res[,3])

############### Figure 4 ################
cv.plot <- ggplot(cv.df,aes(x=mod,y=ll)) + 
  geom_boxplot() +
  theme_bw(base_size=19) + 
  xlab("Model") +
  ylab("Mean log likelihood in test set")

pdf("figures_and_tables/figure4.pdf",width=12,height=5)
cv.plot
dev.off()

############################# Permutation results #############################
load("data_example_results/data_example_full.Rdata")
beta_est_true <- beta_est
group_est_true <- as.integer(as.factor(beta_est_true))

nperm <- 10
perm.est <- matrix(nrow=pL,ncol=nperm)
perm.group <- matrix(nrow=pL,ncol=nperm)
perm.GRI <- matrix(nrow=max(group_est_true),ncol=nperm)

for(i in 1:(nperm)){
  load(paste0("data_example_results/data_example_full_perm",i,".Rdata"))
  perm.est[,i] <- beta_est
  perm.group[,i] <- as.integer(as.factor(beta_est))
  perm.GRI[,i] <- GRI_fun(group_est_true,perm.group[,i])
  #plot(beta_est,beta_est_true)
}
GRI.mean <- rowMeans(perm.GRI)
cat("\n\nGRI_g for each original group, averaged across 10 permutations:\n\n")
cbind(g=1:8,GRI_g=GRI.mean)

dat <- as.data.frame(cbind(perm.est,perm.group))
est.names <- sapply(1:10,function(i) paste0("est.",i))
names(dat)[1:10] <- est.names
group.names <- sapply(1:10,function(i) paste0("group.",i))
names(dat)[11:20] <- group.names
dat$est.orig <- beta_est_true
dat$group.orig <- group_est_true
dat$est.0 <- dat$est.orig
dat$group.0 <- dat$group.orig
dat$Y <- Y
dat.df <- reshape(dat,direction="long",
                  varying=list(c("est.0",est.names),
                               c("group.0",group.names)))
names(dat.df)[4:7] <- c("perm","est","group","v")
dat.df$group.orig <- as.factor(dat.df$group.orig)
dat.df$perm <- dat.df$perm-1
dat.df$perm <- as.factor(dat.df$perm)
dat.df$group <- as.factor(dat.df$group)
change.group <- apply(X=cbind(as.character(dat.df$perm),
                              as.character(dat.df$group.orig),
                              as.character(dat.df$group)),
                      MARGIN=1,FUN=glue_collapse,sep="")
dat.df$change.group <- change.group
n.outcomes <- tapply(change.group,as.factor(change.group),length)
n.cases <- tapply(dat.df$Y,as.factor(change.group),sum)
n.df <- data.frame(change.group=names(n.outcomes),
                   n.outcomes=n.outcomes,
                   n.cases=n.cases)
dat.df <- merge(dat.df,n.df,by="change.group",all.x=T,all.y=F)
dat.df$est.group.orig <- exp(dat.df$est.orig*10)
dat.df$est.group.orig <- as.factor(sprintf("%.3f",dat.df$est.group.orig))
dat.df$est.group <- exp(dat.df$est*10)
dat.df$est.group <- as.factor(sprintf("%.3f",dat.df$est.group))
perm.GRI <- cbind(rep(1,nrow(perm.GRI)),perm.GRI)

##### OR for original groups vs permuted groups, showing number of outcomes
plot.list <- list()
# plotly.list <- list()
plotly.df <- data.frame(matrix(nrow=0,ncol=13))
top.pts <- 10
left.pts <- 10
for(i in 0:10){
  # groups plot
  dat.i <- dat.df[dat.df$perm==i,]
  dat.i$est.group2 <- as.factor(dat.i$est.group)
  dat.i$node_icd9 <- as.character(icd_short_to_decimal(V(tree)$name[dat.i$v+(p-pL)]))
  nodeslist <- tapply(dat.i$node_icd9,dat.i$change.group,
                      glue_collapse,sep=", ",width=80)
  nodeslist <- data.frame(change.group=names(nodeslist),nodes=nodeslist)
  dat.i$v <- NULL
  dat.i$node_icd9 <- NULL
  dat.i$Y <- NULL
  dat.i <- dat.i[!duplicated(dat.i),]
  dat.i <- merge(dat.i,nodeslist,by="change.group")
  if(i > 0){
    plot_title <- paste0("Permutation ",i)
  } else {
    plot_title <- "No permutation" 
  }
  plot.groups <- ggplot(dat.i,aes(x=as.factor(est.group),y=est.group.orig,
                                  label=n.outcomes,text=nodes)) + 
    geom_point(color="white",size=6) +
    # geom_label(size=4,label.size=0,label.padding=unit(0,"lines")) +
    geom_text(size=4) +
    theme_bw() +
    theme(legend.position="none",
          plot.margin = margin(t=top.pts,r=0,b=10,l=left.pts, unit = "pt"),
          plot.title = element_text(hjust = 0.5),
          axis.text=element_text(size=7)) +
    ggtitle(plot_title) +
    xlab("Permuted OR") +
    ylab("Original OR")
  
  if(i > 0){
    # plotly.list[[i]] <- ggplotly(plot.groups,tooltip="text")
    dat.i$perm <- paste0("Permutation ",i)
    names(plotly.df) <- names(dat.i)
    plotly.df <- rbind(plotly.df,dat.i)
  }

  # GRI_g plot
  GRI.df <- data.frame(OR=as.factor(unique(dat.df$est.group.orig)),smc=perm.GRI[,i+1],perm="GRI")
  GRI.df$label <- sprintf("%.3f",GRI.df$smc)
  if(i==0) GRI.df$label <- rep("1.00",nrow(GRI.df))
  plot.GRI <- ggplot(GRI.df,aes(y=OR,x=1)) + 
    geom_tile(aes(fill=GRI),colour="black",size=0.1) +
    geom_label(aes(x=1,y=OR,label=label),size=2,alpha=0.6,label.size=0,label.padding=unit(0.1,"lines")) +
    # geom_text(aes(x=1,y=OR,label=label)) +
    theme_bw() +
    theme(legend.position="none",
          #axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.margin = margin(t=top.pts,r=0,b=0, unit = "pt"),
          plot.title = element_text(hjust = 0.5),
          panel.border=element_blank()) +
    ggtitle(bquote(GRI[g])) +
    xlab("") + 
    scale_fill_gradientn(colours=brewer.pal(max(as.numeric(dat.df$group.orig)),
                                              name="YlGnBu"),
                         limits=c(0,1)) +
    scale_x_continuous(breaks=0.5,labels="")

  # Save plots as pdf
  if(i > 0){
    plot.list[[i]] <- plot.groups + plot.smc +
      plot_layout(ncol=2,widths=c(7,1))
  } 
  # else {
  #   pdf(file = paste0("./figures_and_tables/figureA4_",i,".pdf"),
  #       width=4,height=3.5)
  #   print(plot.groups + plot.smc +
  #           plot_layout(ncol=2,widths=c(7,1)))
  #   dev.off()
  # }
}

## interactive plot
plotly.df$perm <- factor(plotly.df$perm,levels=unique(plotly.df$perm))
plotly.groups <- ggplot(plotly.df,aes(x=as.factor(est.group),y=est.group.orig,
                 label=n.outcomes,text=nodes)) + 
  # geom_point(color="white",size=6) +
  geom_text(size=4) +
  theme_bw() +
  theme(legend.position="none",
        plot.margin = margin(t=top.pts,r=0,b=10,l=left.pts, unit = "pt"),
        axis.text=element_text(size=7)) +
  xlab("Permuted OR") +
  ylab("Original OR") +
  facet_wrap(.~perm,ncol=4,scales="free_x")

# run this to examine which outcomes fall into which groups
ggplotly.groups <- ggplotly(plotly.groups,tooltip="text")
setwd("./figures_and_tables/")
htmlwidgets::saveWidget(as_widget(ggplotly.groups), "permutation_groups.html")
setwd("../")

# wrapping parameters
ncol <- 3
nrow <- ceiling(nperm / ncol)
nblanks <- ncol*nrow - nperm
# add blank plots for wrapping purposes
for(i in 1:nblanks){
  plot.list[[nperm+i]] <- ggplot() + theme_void()
}

# Save as pdf
plot.widths <- 4
plot.heights <- 4
pdf(file = paste0("./figures_and_tables/figureA4.pdf"),
    width=ncol*plot.widths,
    height=nrow*plot.heights)
wrap_plots(plot.list,ncol=3,widths=rep(1,3),heights=rep(1,2))
dev.off()

############################# Sensitivity analysis results #############################
load("data_example_results/data_example_full.Rdata")
beta_est_orig <- exp(final_ss$moretrees_est*10)
beta_group_orig <- as.integer(as.factor(beta_est_orig))
load("data_example_results/data_example_full_sensitivity.Rdata")
beta_est_sens <- exp(final_ss$moretrees_est*10)
beta_group_sens <- as.integer(as.factor(beta_est_sens))
# Load tree
load("simulation_inputs/inputs.Rdata")

## data for plotting
change.group <- apply(X=cbind(as.character(beta_group_orig),
                              as.character(beta_group_sens)),
                      MARGIN=1,FUN=glue_collapse,sep="")
counts <- table(change.group)
counts <- data.frame(change.group=names(counts),n.outcomes=as.integer(counts))
dat <- data.frame(group_orig=beta_group_orig,
                  est_orig=beta_est_orig,
                  group_sens=beta_group_sens,
                  est_sens=beta_est_sens,
                  change.group=change.group)
dat$node_icd9 <- as.character(icd_short_to_decimal(V(tree)$name[(p-pL+1):p]))
nodeslist <- tapply(dat$node_icd9,dat$change.group,
                    glue_collapse,sep=", ",width=40)
nodeslist <- data.frame(change.group=names(nodeslist),nodes=nodeslist)
dat$node_icd9 <- NULL

# for group plot
dat.plt<- dat[!duplicated(dat),]
dat.plt <- merge(dat.plt,counts,by="change.group")
dat.plt <- merge(dat.plt,nodeslist,by="change.group")
dat.plt$lab_orig <- sprintf("%.3f",dat.plt$est_orig)
dat.plt$lab_orig <- factor(dat.plt$lab_orig,levels=sort(unique(dat.plt$lab_orig)))
dat.plt$lab_sens <- sprintf("%.3f",dat.plt$est_sens)
dat.plt$lab_sens <- factor(dat.plt$lab_sens,levels=sort(unique(dat.plt$lab_sens)))

#### Comparison plot ####
top.pts <- 10
left.pts <- 10
# groups plot
plot.groups <- ggplot(dat.plt,aes(x=lab_sens,y=lab_orig,
                                label=n.outcomes,text=nodes)) + 
  geom_point(color="white",size=6) +
  # geom_label(size=4,label.size=0,label.padding=unit(0,"lines")) +
  geom_text(size=4) +
  theme_bw() +
  theme(legend.position="none",
        plot.margin = margin(t=top.pts,r=0,b=10,l=left.pts, unit = "pt"),
        plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=7)) +
  xlab("Sensitivity analysis OR") +
  ylab("Original analysis OR")

# simple matching coefficient plot
smc.sens <- GRI_fun(x=dat$group_orig,y=dat$group_sens)
smc.df <- data.frame(group_orig=1:length(smc.sens),smc=smc.sens)
dat2 <- dat[,c("group_orig","est_orig")]
dat2 <- dat2[!duplicated(dat2),]
smc.df <- merge(smc.df,dat2,by="group_orig")
smc.df$label <- sprintf("%.3f",smc.df$smc)
smc.df$est_lab <- sprintf("%.3f",smc.df$est_orig)
plot.smc <- ggplot(smc.df,aes(y=est_lab,x=1)) + 
  geom_tile(aes(fill=smc),colour="black",size=0.1) +
  geom_label(aes(x=1,y=est_lab,label=label),size=2,alpha=0.6,label.size=0,label.padding=unit(0.1,"lines")) +
  # geom_text(aes(x=1,y=est_lab,label=label)) +
  theme_bw() +
  theme(legend.position="none",
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = margin(t=top.pts,r=0,b=0, unit = "pt"),
        plot.title = element_text(hjust = 0.5),
        panel.border=element_blank()) +
  ggtitle(bquote(GRI[g])) +
  xlab("") + 
  scale_fill_gradientn(colours=brewer.pal(max(as.numeric(smc.df$group_orig)),
                                          name="YlGnBu"),
                       limits=c(0,1)) +
  scale_x_continuous(breaks=0.5,labels="")

### save as pdf
pdf(file="./figures_and_tables/figureA5.pdf",width=4.2,height=4)
plot.groups + plot.smc + plot_layout(ncol=2,widths=c(7,1))
dev.off()

### interactive plot ###
plotly.groups <- ggplot(dat.plt,aes(y=lab_orig,x=lab_sens,
                                      label=n.outcomes,text=nodes)) + 
  geom_text(size=4) +
  theme_bw() +
  theme(legend.position="none",
        axis.text=element_text(size=10)) +
  xlab("Sensitivity analysis OR") +
  ylab("Original analysis OR") 

# run this to examine which outcomes fall into which groups
ggplotly.groups <- ggplotly(plotly.groups,tooltip="text",width=500,height=500) %>% 
  layout(margin=list(l = -10), yaxis=list(tickprefix="    "))
setwd("./figures_and_tables/")
htmlwidgets::saveWidget(as_widget(ggplotly.groups), "sensitivity_analysis_groups.html")
setwd("../")

############################### Full data results ###############################
load("data_example_results/data_example_full.Rdata")

### Compute grouped estimates and CIs ###
beta_est <- final_ss$moretrees_est
VI_params <- final_ss$VI_params
beta_groups <- groups
tree <- groups.calc.fun(tree,beta_groups,beta_est,VI_params)
beta_indiv_est <- indiv.beta.ci.calc(VI_params,ancestors,pL,p)
V(tree)$beta_indiv <- numeric(length=p)
V(tree)$beta_indiv[V(tree)$leaf] <- beta_indiv_est$beta_indiv
V(tree)$beta_indiv_cil <- numeric(length=p)
V(tree)$beta_indiv_cil[V(tree)$leaf] <- beta_indiv_est$cil_indiv
V(tree)$beta_indiv_ciu <- numeric(length=p)
V(tree)$beta_indiv_ciu[V(tree)$leaf] <- beta_indiv_est$ciu_indiv

############## Table 1 accompanying figure ################
mult <- 10
beta_est_mult <- beta_est*mult
l.big <- layout.reingold.tilford(tree, root=root.node)
l.big[,2] <- -40*l.big[,2]
l.big[,1] <- 4*l.big[,1]
l.big[488:499,1] <- l.big[488:499,1] + 6
l.big[149:151,1] <- l.big[149:151,1] + c(8,16,24)
beta_out <- final_ss$moretrees_est
cols <- rainbow(7, alpha = 1)
col1 <- rep(NA,p[1])
col1[V(tree)$leaf] <- groups
plot.groups <- list()
for(i in 1:6){
  plot.groups[[i]] <- which(col1==i)
  col1[col1==i] <- cols[i]
}
plot.groups[[7]] <- which(col1==8)
col1[col1==8] <- cols[7]
col1[col1==7] <- "grey"
col1[is.na(col1)] <- "white"

pdf("./figures_and_tables/table1_tree.pdf",width=10,height=4)
plot(tree, layout = l.big, margin=c(1,1,1,1), rescale=F, xlim=c(min(l.big[,1]),max(l.big[,1])),
     ylim=c(max(l.big[,2]),min(l.big[,2])),
     vertex.label=NA,
     edge.arrow.size=0,margin=c(0,0,0,0),vertex.size=300,
     vertex.color=col1,edge.width=0.5,vertex.frame.color=NA,
     mark.groups=plot.groups,mark.col=rainbow(7, alpha = 0.4),
     mark.border = rainbow(length(plot.groups), alpha = 0.8), mark.expand = 1000)
dev.off()

cols2 <- c(cols[1:6],"grey",cols[7])
cols.box <- rainbow(7, alpha = 0.4)
cols.box <- c(cols.box[1:6],NA,cols.box[7])
cols.border <- rainbow(7, alpha = 0.8)
cols.border <- c(cols.border[1:6],NA,cols.border[7])
group.names <- c("Group 1","Group 2","Group 3","Group 4","Group 5","Group 6","Group 7","Group 8")
pdf("figures_and_tables/table1_legend.pdf",width=14,height=3)
plot(c(-1,1), c(-1,1), type = "n", ann = F, axes = F)
legend('center',legend=group.names,col=cols.box,pch=19,pt.cex=3,y.intersp=2,x.intersp=1.5,bty="n",horiz=T,text.width=rep(0.15,8))
legend('center',legend=group.names,col=cols.border,pch=1,pt.cex=3,pt.lwd=2,y.intersp=2,x.intersp=1.5,bty="n",horiz=T,text.width=rep(0.15,8))
legend('center',legend=group.names,col=cols2,pch=19,y.intersp=2,x.intersp=1.5,bty="n",horiz=T,text.width=rep(0.15,8))
dev.off()

col2 <- rep("black",p)
col2[V(tree)$leaf] <- "red"
l.big2 <- l.big
l.big2[488:499,1] <- l.big2[488:499,1] + 6
l.big2[,2] <- l.big2[,2]*2
pdf("figures_and_tables/figure2.pdf",width=10,height=4)
plot(tree, layout = l.big, margin=c(1,1,1,1), rescale=F, xlim=c(min(l.big[,1]),max(l.big[,1])),
     ylim=c(max(l.big[,2]),min(l.big[,2])),
     vertex.label=NA,edge.arrow.size=0.05,margin=c(0,0,0,0),vertex.size=300,
     vertex.color=col2,edge.width=0.5,vertex.frame.color=NA)
dev.off()

################# Interactive Trees ##########################

#### Prune the tree so that all siblings with same values are merged ####

mult <- 10 # show results in units of 10 micrograms per cubic metre
digits <- 3 # number of digits to display

{
  # Step 1: collapse to parent node when all children have same beta value
  trPrune <- tree
  V(trPrune)$outcomes <- sapply(V(trPrune)$name,list)
  V(trPrune)$nsamp <- numeric(length=length(V(trPrune))) + 0
  V(trPrune)$nsamp[V(trPrune)$leaf] <- Y
  nodes <- V(trPrune)
  iter <- 1
  order <- 1000
  while(iter < length(nodes)){
    v <- nodes[iter]
    children <- as.numeric(ego(trPrune,mode="out",nodes=v,order=order,mindist=1)[[1]])
    beta.children <- V(trPrune)$beta_grouped[children]
    beta.children.cil <- V(trPrune)$beta_grouped_cil[children]
    beta.children.ciu <- V(trPrune)$beta_grouped_ciu[children]
    if(sum(is.na(beta.children))==0 & length(unique(beta.children))==1){
      # If all children have same beta value, merge them into parent
      V(trPrune)$outcomes[v][[1]] <- unique(c(unlist(vertex.attributes(trPrune,children)$outcomes),names(v)))
      V(trPrune)$leaf[v] <- TRUE
      V(trPrune)$beta_grouped[v] <- unique(beta.children)
      V(trPrune)$beta_grouped_ciu[v] <- unique(beta.children.ciu)
      V(trPrune)$beta_grouped_cil[v] <- unique(beta.children.cil)
      V(trPrune)$nsamp[v] <- sum(V(trPrune)$nsamp[children])
      trPrune <- delete.vertices(trPrune,children)
      iter <- 1
    } else {
      iter <- iter + 1
    }
    nodes <- V(trPrune)
  }
  
  # Step 2: merge leaf siblings with same values
  iter <- 1
  leaves <- as.numeric(V(trPrune)[V(trPrune)$leaf])
  V(trPrune)$parent <- numeric(length(V(trPrune))) + NA # NA means no parent (root node)
  for(v in V(trPrune)){
    children <- ego(trPrune,order=1,mindist=1,nodes=v,mode="out")[[1]]
    if(length(children)>0){
      V(trPrune)$parent[children] <- as.numeric(v)
    }
  }
  leaves <- V(trPrune)$name[V(trPrune)$leaf]
  sibs.leaves <- vertex.attributes(trPrune,leaves)$parent
  beta.leaves <- vertex.attributes(trPrune,leaves)$beta_grouped
  sibscols <- cbind(as.character(sibs.leaves),as.character(beta.leaves))
  groupings <- as.numeric(as.factor(apply(sibscols,1,glue_collapse,sep=",")))
}
# groupings is a vector that uniquely identifies groups of leaf siblings with the same beta value.
# These need to be merged.
while(length(leaves)>0){
  gr <- groupings[1]
  nodes <- leaves[groupings==gr]
  if(length(nodes)>1){
    newOutcomes <- unlist(vertex.attributes(trPrune,nodes)$outcomes)
    newOutcomes <- unique(c(newOutcomes,nodes))
    n_new <- sum(V(trPrune)$nsamp[names(V(trPrune)) %in% nodes])
    vertex_attr(trPrune,name="outcomes",index=nodes[1])[[1]] <- newOutcomes
    trPrune <- delete.vertices(trPrune,nodes[-1])
    vertex_attr(trPrune,name="nsamp",index=nodes[1]) <- n_new
  }
  leaves <- leaves[groupings!=gr]
  groupings <- groupings[groupings!=gr]
}

edgelist2 <- as.data.frame(get.edgelist(trPrune),stringsAsFactors=F)
names(edgelist2) <- c("parent","child")
edgelist2 <- rbind(c(NA,root.node),edgelist2)
edgelist2$beta <- get.vertex.attribute(trPrune,"beta_grouped",edgelist2$child)*mult
edgelist2$beta_cil <- get.vertex.attribute(trPrune,"beta_grouped_cil",edgelist2$child)*mult
edgelist2$beta_ciu <- get.vertex.attribute(trPrune,"beta_grouped_ciu",edgelist2$child)*mult
edgelist2$nsamp <- get.vertex.attribute(trPrune,"nsamp",edgelist2$child)
maxbeta <- max(abs(final_ss$moretrees_est))
cols2 <- colorRamp2(c(-maxbeta, 0, maxbeta), c("red", "white", "blue"))
colsfun <- function(beta){
  if(is.na(beta)){
    return("black")
  } else {
    return(cols2(as.numeric(beta)))
  }
} 
edgelist2$col <- sapply(edgelist2$beta,colsfun)
edgelist2$outcomes <- get.vertex.attribute(trPrune,"outcomes",edgelist2$child)
leaves <- get.vertex.attribute(trPrune,"leaf",edgelist2$child)
edgelist2$leaf <- leaves
edgelist2$col[!leaves] <- "green"
edgelist2$nodesize <- 1
edgelist2$nodesize[!leaves] <- 0.001

# html tooltip
source("scripts/collapsibleTreeNetwork_modified.R")
edgelist2$tooltip <- character(length=nrow(edgelist2))
for(i in 1:nrow(edgelist2)){
  edgelist2$tooltip[i] <- expand_html(edgelist2[i,],digits=digits,tr=tree)
}
edgelist2_parent <- as.character(icd_short_to_decimal(edgelist2$parent))
edgelist2$parent[!is.na(edgelist2_parent)] <- edgelist2_parent[!is.na(edgelist2_parent)]
edgelist2_child <- as.character(icd_short_to_decimal(edgelist2$child))
edgelist2$child[!is.na(edgelist2_child)] <- edgelist2_child[!is.na(edgelist2_child)]

# Plot
pl <-   collapsibleTreeNetwork2(edgelist2,attribute="beta",aggFun=identity, fill="col",
                               tooltipHtml="tooltip",nodeSize="nodesize",nodeSizeAggFun=identity,
                               nodeSizeScaleFun=max,width=800,height=1000,fontSize=14)
setwd("figures_and_tables/")
htmlwidgets::saveWidget(as_widget(pl), "data_example_collapsed_tree.html")
setwd("../")

#### Individual level tree ####
edgelist3 <- as.data.frame(get.edgelist(tree),stringsAsFactors=F)
names(edgelist3) <- c("parent","child")
edgelist3 <- rbind(c(NA,root.node),edgelist3)
edgelist3$beta <- get.vertex.attribute(tree,"beta_indiv",edgelist3$child)*mult
edgelist3$beta_cil <- get.vertex.attribute(tree,"beta_indiv_cil",edgelist3$child)*mult
edgelist3$beta_ciu <- get.vertex.attribute(tree,"beta_indiv_ciu",edgelist3$child)*mult
V(tree)$nsamp <- numeric(p)
V(tree)$nsamp[V(tree)$leaf] <- Y
edgelist3$nsamp <- get.vertex.attribute(tree,"nsamp",edgelist3$child)
maxbeta <- max(abs(V(tree)$beta_indiv))
cols3 <- colorRamp2(c(-maxbeta, 0, maxbeta), c("red", "white", "blue"))
edgelist3$col <- sapply(edgelist3$beta,colsfun)
V(tree)$outcomes <- sapply(V(tree)$name,list)
edgelist3$outcomes <- get.vertex.attribute(tree,"outcomes",edgelist3$child)
leaves <- get.vertex.attribute(tree,"leaf",edgelist3$child)
edgelist3$leaf <- leaves
edgelist3$col[!leaves] <- "green"
edgelist3$nodesize <- 1
edgelist3$nodesize[!leaves] <- 0.001
edgelist3$nodesize <- edgelist3$nodesize

# html tooltip
edgelist3$tooltip <- character(length=nrow(edgelist3))
for(i in 1:nrow(edgelist3)){
  edgelist3$tooltip[i] <- expand_html(edgelist3[i,],digits=digits,tr=tree)
}
edgelist3_parent <- as.character(icd_short_to_decimal(edgelist3$parent))
edgelist3$parent[!is.na(edgelist3_parent)] <- edgelist3_parent[!is.na(edgelist3_parent)]
edgelist3_child <- as.character(icd_short_to_decimal(edgelist3$child))
edgelist3$child[!is.na(edgelist3_child)] <- edgelist3_child[!is.na(edgelist3_child)]

# Plot
silly_function <- function(x) 1
pl2 <-   collapsibleTreeNetwork2(edgelist3,attribute="beta",aggFun=identity, fill="col",
                                tooltipHtml="tooltip",nodeSize="nodesize",nodeSizeAggFun=identity,
                                nodeSizeScaleFun=silly_function,width=800,height=1000,fontSize=14)
#htmlwidgets::saveWidget(as_widget(pl2), "/Users/emt380/Documents/emgthomas.github.io/moretrees_plots/data_example_individual_tree.html")
setwd("./figures_and_tables/")
htmlwidgets::saveWidget(as_widget(pl2), "data_example_individual_tree.html")
setwd("../")

################### Table 1 #######################

latex_groups <- character(length=max(beta_groups))
small_table <- data.frame(group=character(max(beta_groups)),
                          codes=character(max(beta_groups)),
                          num_cases=character(max(beta_groups)),
                          smc=character(max(beta_groups)),
                          OR_moretrees=character(max(beta_groups)),
                          OR_ml=character(max(beta_groups)),
                          stringsAsFactors = F)
compare_models <- matrix(nrow=max(beta_groups),ncol=6)
for(g in 1:max(beta_groups)){
  ancestors_g <- unique(unlist(ancestors[which(beta_groups==g) + (p-pL)]))
  tree_g <- induced_subgraph(tree,V(tree)[ancestors_g])
  nodes <- names(V(tree_g))
  leaves <- names(V(tree_g))[V(tree_g)$leaf]
  beta_moretrees <- V(tree_g)$beta_grouped[V(tree_g)$leaf][1]*mult
  ci_moretrees <- c(V(tree_g)$beta_grouped_cil[V(tree_g)$leaf][1],V(tree_g)$beta_grouped_ciu[V(tree_g)$leaf][1])*mult
  beta_ml <- beta.ml.groups[g,1]*mult
  ci_ml <- beta.ml.groups[g,2:3]*mult
  nsamp <- Y[beta_groups==g]
  compare_models[g,] <- c(beta_moretrees,beta_ml,ci_moretrees,as.numeric(ci_ml))
  # big table
  latex_groups[g] <- expand_groups_latex(g,nodes,beta_moretrees,ci_moretrees,nsamp,digits=digits,tree_g,indent.multiplier = 10)
  # small table
  frmt <- paste0("%.",digits,"f")
  beta_moretrees_frmt <- sprintf(frmt,exp(beta_moretrees))
  beta_moretrees_ciu_frmt <- sprintf(frmt,exp(ci_moretrees[2]))
  beta_moretrees_cil_frmt <- sprintf(frmt,exp(ci_moretrees[1]))
  beta_ml_frmt <- sprintf(frmt,exp(beta_ml))
  beta_ml_ciu_frmt <- sprintf(frmt,exp(ci_ml[2]))
  beta_ml_cil_frmt <- sprintf(frmt,exp(ci_ml[1]))
  smc_frmt <- sprintf(frmt,smc.mean[g])
  if(g != 7){
    codes_g <- as.character(icd_short_to_decimal(leaves))
  } else {
    codes_g <- "All 420 remaining codes"
  }
  small_table[g,] <- c(g,format(paste(codes_g,collapse=", "),big.mark=","),
                       format(sum(nsamp),big.mark=","),
                       smc_frmt,
                       paste0(beta_moretrees_frmt," (",beta_moretrees_cil_frmt,",",beta_moretrees_ciu_frmt,")"),
                       paste0(beta_ml_frmt," (",beta_ml_cil_frmt,",",beta_ml_ciu_frmt,")"))
}

# write big table
write(paste("\\renewcommand\\arraystretch{0.6}\\begin{longtable}{p{\\textwidth}} \\hline",paste(latex_groups,sep="",collapse=" \n \\hline \n \n "),"\\end{longtable}",sep="\n \n"),file="./figures_and_tables/tableC_supplementary_material.tex")

# write small table
row.names(small_table) <- NULL
small_xtable <- xtable(small_table,align=c("l","c","p{6cm}","c","c","c","c"),digits=3,
                       display=c("d","d","s","d","d","f","f"))
names(small_xtable) <- c("Group","ICD9 codes","#Cases","Permuted SMC","OR (95%CI) MOReTreeS","OR (95% CI) Maximum Likelihood")

write(print(small_xtable,floating=FALSE,include.rownames = FALSE),file="./figures_and_tables/table1.tex")

########################## Conditional coverage analysis #############################

# Appending simulation blocks
VI_files <- list.files(path="simulation_results/coverage/",pattern="VI_coverage_block*")
VI_new_file <- "simulation_results/VI_conditional_coverage.csv"
file.create(VI_new_file)
write.table(rbind(c("sim","mu_gamma","sigma2_gamma","u_s")), file = VI_new_file, row.names=FALSE, col.names=FALSE, sep=",",append=TRUE)
i <- 0
for(fn in VI_files){
  sims <- read.csv(file=paste0("simulation_results/coverage/",fn),header=F,row.names=NULL)
  sims[,1] <- sims[,1] + i
  i <- max(sims[,1])
  write.table(sims, file=VI_new_file, row.names =FALSE, col.names = FALSE,sep = ",", append = TRUE)
}

beta_ML_files <- list.files(path="simulation_results/coverage/",pattern="beta_ML_coverage_block*")
beta_ML_new_file <- "simulation_results/beta_ML_conditional_coverage.csv"
file.create(beta_ML_new_file)
write.table(rbind(c("sim","beta_ML","ci_l","ci_u")), file = beta_ML_new_file, row.names=FALSE, col.names=FALSE, sep=",",append=TRUE)
i <- 0
for(fn in beta_ML_files){
  sims <- read.csv(file=paste0("simulation_results/coverage/",fn),header=F,row.names=NULL)
  sims[,1] <- sims[,1] + i
  i <- max(sims[,1])
  write.table(sims, file=beta_ML_new_file, row.names =FALSE, col.names = FALSE,sep = ",", append = TRUE)
}

rm(sims)

####### Get coverage per group #######

# Load simulation results
VIsims <- read.csv("simulation_results/VI_conditional_coverage.csv",header = T)
betaMLsims <- read.csv("simulation_results/beta_ML_conditional_coverage.csv",header = T)
# Load full data results
load("data_example_results/data_example_full.Rdata",verbose = T)
beta_true <- beta_est # used the estimated betas as true betas in this simulation
groups_true <- groups

# Load tree
load("simulation_inputs/inputs.Rdata")

# Get ancestor matrix A
A <- t(as_adj(tree,sparse = T))
A <- expm(A)
A[A>0] <- 1 # row A[v,] indicates which nodes are ancestors of v

# Modify it so that only non-zero nodes are equal to one
nonzero <- exp(loglogit(final_ss$VI_params$u_s)) >= 0.5
A.mod <- A %*% diag(nonzero)
A.conv <- A.mod[(p-pL+1):p,] # matrix for conversion: beta <- A.conv %*% gamma

# check it worked
beta_true <- as.numeric(A.conv %*% final_ss$VI_params$mu_gamma)
all.equal(beta_true,beta_est) # yep

# list to store coverage results
coverage_VI <- rep(list(numeric(0)),8)
coverage_ML <- rep(list(numeric(0)),8)
pbias_VI <- rep(list(numeric(0)),8)
pbias_ML <- rep(list(numeric(0)),8)
# get binary coverage indicator for each simulation
for(i in 1:max(VIsims$sim)){
  VIsims_i <- VIsims[VIsims$sim==i,]
  betaMLsims_i <- betaMLsims[betaMLsims$sim==i,]
  # Get CIs for simulation
  nonzero <- exp(loglogit(VIsims_i$u_s)) >= 0.5
  A.mod <- A %*% diag(nonzero)
  A.conv <- A.mod[(p-pL+1):p,] # matrix for conversion: beta <- A.conv %*% gamma
  beta_est <- as.numeric(A.conv %*% VIsims_i$mu_gamma)
  groups_est <- as.integer(as.factor(beta_est))
  gtab <- table(groups_est,groups_true)
  for(g in 1:max(groups_true)){ # for each true group
    if(sum(gtab[,g]>0)==1){ # did true group get split into multiple estimated groups?
      # If no...
      h <- as.integer(which(gtab[,g]>0)) # which estimated group contains all outcomes in group g?
      if(sum(gtab[h,]>0)==1){ # does this group contain other outcomes, outside group g?
        # If no, groups g and h are identical. Proceed with estimating coverage.
        
        # Did the moretrees estimate cover the true coefficient?
        beta_true_g <- unique(beta_true[groups_true==g]) # what was the true beta for this group?
        v <- which(groups_est==h)[1] # select a representative outcome from the group (results will be same for any such v)
        beta_est_g <- unique(beta_est[v]) # get estimated beta for group h
        sd_est_g <- sqrt(as.numeric(A.conv[v,] %*% VIsims_i$sigma2_gamma)) # get sd for group h
        ci_l <- beta_est_g - qnorm(0.975)*sd_est_g # lower bound of CI
        ci_u <- beta_est_g + qnorm(0.975)*sd_est_g # upper bound of CI
        covered <- (ci_l <= beta_true_g) & (ci_u >= beta_true_g) # does the CI cover the true beta?
        coverage_VI[[g]] <- c(coverage_VI[[g]],covered) # store in list
        pbias_VI[[g]] <- c(pbias_VI[[g]],abs(beta_est_g-beta_true_g)/abs(beta_true_g))
        
        # Did the ML estimate cover the true coefficient?
        beta_ML_g <- unique(betaMLsims_i$beta_ML[groups_true==g])
        ci_ML_l <- unique(betaMLsims_i$ci_l[groups_true==g])
        ci_ML_u <- unique(betaMLsims_i$ci_u[groups_true==g])
        covered_ML <- (ci_ML_l <= beta_true_g) & (ci_ML_u >= beta_true_g)
        coverage_ML[[g]] <- c(coverage_ML[[g]],covered_ML) 
        pbias_ML[[g]] <- c(pbias_ML[[g]],abs(beta_ML_g-beta_true_g)/abs(beta_true_g))
      }
    }
  }
  if(i %% 100 == 0) print(i)
}

# what's our conditional sample size for each group?
nsamp_g <- sapply(coverage_ML,length)

########### Table A.2 #############
# compute group-specific coverage
cov_groups_VI <- sapply(coverage_VI,mean)
cov_groups_ML <- sapply(coverage_ML,mean)
cov_tab <- data.frame("Group"=1:8,"ssMOReTreeS"=cov_groups_VI*100,"Maximum likelihood"=cov_groups_ML*100)
print(xtable(cov_tab),include.rownames=F)

# compute group-specific percentage bias
pbias_groups_VI <- sapply(pbias_VI,mean)
pbias_groups_ML <- sapply(pbias_ML,mean)
cbind(pbias_groups_ML,pbias_groups_VI)

# remove csv files because they are too large for github
file.remove(beta_ML_new_file)
file.remove(VI_new_file)
