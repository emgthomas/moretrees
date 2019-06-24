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

############################### Cross-validation results ##############################

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
plot(tp,ll.tp,type="l")

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
  ngroups[i] <- length(unique(as.numeric(beta.t0)))
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
tp.marks <- c(tp[1],min(tp.max),max(tp.max),tp[length(tp)])
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

############### Figure A3 ################
pdf("figures_and_tables/figureA3_a.pdf",width=8,height=5)
plot_tp
dev.off()

pdf("figures_and_tables/figureA3_b.pdf",width=8,height=5)
plot_ngroups
dev.off()

### Load other CV results ###
nmods <- 9
cv.res <- as.data.frame(matrix(nrow=nfolds,ncol=nmods+1))
names(cv.res) <- c("fold",
                   "ssMOReTreeS\n Collapsed",
                   "ssMOReTreeS\n Collapsed\n Max Delta",
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
  cv.res[i,] <- c(ll.cv[1,1:2],ll.max.folds[1,i],as.numeric(ll.cv[1,3:length(ll.cv)]))
}
require(data.table)
cv.res$sim_groups <- NULL
cv.df <- melt(cv.res,id.vars="fold",measure.vars=2:(nmods),variable.name="mod",value.name="ll")

# Cross validated predictive accuracy
sort(colMeans(cv.res[,2:ncol(cv.res)]),decreasing = T)

# Winner excluding MOReTreeS individual and MOReTreeS collapsed delta max
apply(cv.res[,-c(1,3,4)],1,which.max)
# How often does MOReTreeS individual beat MOReTrees collapsed?
sum(cv.res[,2] <= cv.res[,4])
# How often does MOReTreeS collpsed delta max beat MOReTrees collapsed?
sum(cv.res[,2] <= cv.res[,3])

############### Figure 4 ################
cv.plot <- ggplot(cv.df,aes(x=mod,y=ll)) + 
  geom_boxplot() +
  theme_bw(base_size=19) + 
  xlab("Model") +
  ylab("Mean log likelihood in test set")

pdf("figures_and_tables/figure4.pdf",width=14,height=5)
cv.plot
dev.off()


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
require(RColorBrewer)
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
require(glue)

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
  groupings <- as.numeric(as.factor(apply(sibscols,1,collapse,sep=",")))
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

require(collapsibleTree)
require(circlize)
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
require(plotly)
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
  if(g != 7){
    codes_g <- as.character(icd_short_to_decimal(leaves))
  } else {
    codes_g <- "All 420 remaining codes"
  }
  small_table[g,] <- c(g,format(paste(codes_g,collapse=", "),big.mark=","),format(sum(nsamp),big.mark=","),paste0(beta_moretrees_frmt," (",beta_moretrees_cil_frmt,",",beta_moretrees_ciu_frmt,")"),paste0(beta_ml_frmt," (",beta_ml_cil_frmt,",",beta_ml_ciu_frmt,")"))
}

# write big table
write(paste("\\renewcommand\\arraystretch{0.6}\\begin{longtable}{p{\\textwidth}} \\hline",paste(latex_groups,sep="",collapse=" \n \\hline \n \n "),"\\end{longtable}",sep="\n \n"),file="./figures_and_tables/tableC_supplementary_material.tex")

# write small table
require(xtable)
row.names(small_table) <- NULL
small_xtable <- xtable(small_table,align=c("l","c","p{6cm}","c","c","c"),digits=3,
                       display=c("d","d","s","d","f","f"))
names(small_xtable) <- c("Group","ICD9 codes","#Cases","OR (95%CI) MOReTreeS","OR (95% CI) Maximum Likelihood")

write(print(small_xtable,floating=FALSE,include.rownames = FALSE),file="./figures_and_tables/table1.tex")



