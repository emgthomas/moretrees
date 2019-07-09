indiv.beta.calc <- function(VI_params,ancestors,pL,p){
  
  ## Inputs ##
  
  # VI_params = a list of variational parameters, as returned from VI_binary_ss
  
  ## Outputs ##
  
  # beta_est, a numeric vector of length pL (the number of outcomes) containing individual estimates
  # of the log odds ratio for each outcome
  
  ## Code ##
  
  pi_s <- 1/(1+exp(-VI_params$u_s))
  mu_gamma <- VI_params$mu_gamma
  beta_est <- numeric(length=pL)
  for(v in 1:pL){
    u <- ancestors[[v+p-pL]]
    beta_est[v] <- sum(mu_gamma[u]*pi_s[u])
  }
  return(beta_est)
}

indiv.beta.sd.calc <- function(idx,VIsims,ancestors,pL,p){
  pi_s <- 1/(1+exp(-VIsims$u_s[idx]))
  mu_gamma <- VIsims$mu_gamma[idx]
  sigma2_gamma <- VIsims$sigma2_gamma[idx]
  beta_est <- numeric(length=pL)
  sd_est <- numeric(length=pL)
  for(v in 1:pL){
    u <- ancestors[[v+p-pL]]
    beta_est[v] <- sum(mu_gamma[u]*pi_s[u])
    sigma2_gamma_u <- sigma2_gamma[u]
    mu_gamma_u <- mu_gamma[u]
    pi_s_u <- pi_s[u]
    sd_est[v] <- sqrt(sum(pi_s_u*(sigma2_gamma_u + (1-pi_s_u)*mu_gamma_u^2)))
  }
  return(data.frame(beta_indiv=beta_est,sd_indiv=sd_est,sim=VIsims$sim[idx[1:pL]]))
}

beta.var.calc <- function(mu_gamma,sigma2_gamma,u_s,ancestors,pL,p){
  pi_s <- 1/(1+exp(u_s))
  # beta_est <- numeric(length=pL)
  var_est <- numeric(length=pL)
  for(v in 1:pL){
    u <- ancestors[[v+p-pL]]
    sigma2_gamma_u <- sigma2_gamma[u]
    mu_gamma_u <- mu_gamma[u]
    pi_s_u <- pi_s[u]
    var_est[v] <- sum(pi_s_u*(sigma2_gamma_u + (1-pi_s_u)*mu_gamma_u^2))
  }
  return(var_est)
}

gamma.sim.fun <- function(u,mu_gamma,sigma2_gamma,n.sim) rnorm(n.sim,mu_gamma[u],sigma2_gamma[u])

indiv.beta.ci.calc <- function(VI_params,ancestors,pL,p,n.sim=1000){
  pi_s <- 1/(1+exp(-VI_params$u_s))
  mu_gamma <- VI_params$mu_gamma
  sigma2_gamma <- VI_params$sigma2_gamma
  ci_l <- numeric(length=pL)
  ci_u <- numeric(length=pL)
  for(v in 1:pL){
    u <- ancestors[[v+p-pL]]
    # Estimate
    beta_est[v] <- sum(mu_gamma[u]*pi_s[u])
    # Simulate credible interval
    s_sim <- sapply(pi_s[u],rbinom,size=1,n=n.sim)
    gamma_sim <- sapply(u,gamma.sim.fun,mu_gamma=mu_gamma,sigma2_gamma=sigma2_gamma,n.sim=n.sim)
    beta_sim <- rowSums(s_sim*gamma_sim)
    ci_l[v] <- quantile(beta_sim,0.025)
    ci_u[v] <- quantile(beta_sim,0.975)
  }
  return(data.frame(beta_indiv=beta_est,cil_indiv=ci_l,ciu_indiv=ci_u))
}

groups.calc.fun <- function(tree,beta_groups,beta_est,VI_params){
  s <- 1/(1+exp(-VI_params$u_s)) >= 0.5
  V(tree)$beta_grouped <- numeric(length(V(tree)))
  V(tree)$beta_grouped_cil <- numeric(length(V(tree)))
  V(tree)$beta_grouped_ciu <- numeric(length(V(tree)))
  for(g in unique(beta_groups)){
    # Find common ancestors
    common.ancestors <- Reduce(intersect,ancestors[which(beta_groups==g)+p-pL])
    active.ancestors <- common.ancestors[s[common.ancestors]]
    # Estimate beta for the group
    beta_est_g <- sum(VI_params$mu_gamma[active.ancestors])
    # Get standard error
    sd_est_g <- sqrt(sum(VI_params$sigma2_gamma[active.ancestors]))
    # Put results into tree
    V(tree)$beta_grouped[which(beta_groups==g)+p-pL] <- beta_est_g
    V(tree)$beta_grouped_cil[which(beta_groups==g)+p-pL] <- beta_est_g + qnorm(p=0.025)*sd_est_g
    V(tree)$beta_grouped_ciu[which(beta_groups==g)+p-pL] <- beta_est_g + qnorm(p=0.975)*sd_est_g
  }
  return(tree)
}

explainer.latex <- function(df){
  if(df["leaf"]==1){
    nsamp_frmt <- as.character(format(as.numeric(df["nsamp"]),big.mark=","))
    glue_collapse(c("\\-\\ \\hspace{",df["indent"],"pt}\\footnotesize{-- {\\color{ForestGreen} \\textbf{",df["icd9_decimal"],"}}: ",df["explainer"]," (n=",nsamp_frmt,")} \\\\ "))
  } else {
    glue_collapse(c("\\-\\ \\hspace{",df["indent"],"pt}\\footnotesize{-- ",df["icd9_decimal"],": ",df["explainer"],"} \\\\ "))
  }
}

expand_groups_latex <- function(g,nodes,beta,ci,nsamp,digits=3,tree_g,indent.multiplier=10){
  # leaf_nodes is a list of leaf nodes included in this group
  # beta is the estimated log odds ratio
  # ci is the corresponding credible interval
  # nsamp is the number of 
  
  frmt <- paste0("%.",digits,"f")
  beta_frmt <- sprintf(frmt,exp(beta))
  beta_ciu_frmt <- sprintf(frmt,exp(ci[2]))
  beta_cil_frmt <- sprintf(frmt,exp(ci[1]))
  nsamp_total <- format(sum(nsamp),big.mark=",")
  chapslist <- sapply(icd9_chapters,function(chap) paste(chap[1],"to",chap[2]))
  chapsnames <- names(icd9_chapters)
  subchapslist <- sapply(icd9_sub_chapters,function(chap) paste(chap[1],"to",chap[2]))
  subchapsnames <- names(icd9_sub_chapters)
  
  #codes <- unlist(V(tr)$outcomes[leaf_nodes+(p-pL)])
  which.chaps <- nodes %in% chapslist
  which.subchaps <- nodes %in% subchapslist
  which.codes <- !(which.chaps | which.subchaps)
  if(sum(which.chaps | which.subchaps | which.codes) != length(nodes)){
    stop("stuffed up somewhere")
  }
  explainer <- character(length=length(nodes))
  icd9_decimal <- character(length=length(nodes))
  if(sum(which.chaps > 0)){
    explainer[which.chaps] <- sapply(nodes[which.chaps], function(code) chapsnames[chapslist == code])
    icd9_decimal[which.chaps] <- nodes[which.chaps]
  }
  if(sum(which.subchaps > 0)){
    explainer[which.subchaps] <- sapply(nodes[which.subchaps], function(code) subchapsnames[subchapslist == code])
    icd9_decimal[which.subchaps] <- nodes[which.subchaps]
  }
  explainer[which.codes] <- sapply(nodes[which.codes],icd_explain,condense=F,brief=F)
  icd9_decimal[which.codes] <- as.character(icd_short_to_decimal(nodes[which.codes]))
  
  # Now get the shape of the 'tree' for this set of nodes
  comp <- components(tree_g)$membership
  ncodes <- sum(V(tree_g)$leaf)
  latex_out <- character()
  for(co in unique(comp)){
    # Get subtree formed by a connected component
    subtr.comp <- induced_subgraph(tree_g,V(tree_g)[comp == co])
    n.comp <- length(V(subtr.comp))
    if(n.comp > 1){
      # Find root node
      root.subtr <- V(subtr.comp)[ego_size(subtr.comp,mode="in",mindist=1,order=1)==0]
      # ordering of nodes using data.tree package
      subtr.comp.data <- data.tree:::FromDataFrameNetwork(as.data.frame(get.edgelist(subtr.comp)))
      codes.order <- subtr.comp.data$Get("name")
      codes.extract <- sapply(codes.order,function(code) which(nodes==code))
      # Amount of indentation (based on distance from root)
      indent <- as.numeric(distances(subtr.comp,v=root.subtr,to=codes.order,mode="out"))*indent.multiplier
      # Put it all in a data frame
      paste.df <- data.frame(icd9_decimal=icd9_decimal[codes.extract],explainer=explainer[codes.extract],
                             indent=indent,leaf=as.numeric(V(subtr.comp)$leaf[codes.extract]),nsamp=V(subtr.comp)$nsamp[codes.extract],row.names=NULL,stringsAsFactors = FALSE)
      # collapse
      latex_out <- c(latex_out,glue_collapse(apply(paste.df,1,explainer.latex)))
    } else {
      latex_out <- c(latex_out,paste0("--",icd9_decimal[comp==co],": ",explainer[comp==co],"\\\\"))
    }
  }
  latex_out <- paste0("\\textbf{\\emph{Group ",g,"}}\\\\ \n\\textbf{Odds Ratio (95\\% Credible Interval) = ",beta_frmt," (",beta_cil_frmt,",",beta_ciu_frmt,")} \\\\ \\textbf{Number of Cases = ",nsamp_total,"} \\\\ ","\\textbf{Group contains the following ",ncodes," ICD9 codes:} \\\\ ",glue_collapse(latex_out))
  return(latex_out)
}

rmse.fun <- function(beta.est,beta) sqrt(mean((beta-beta.est)^2,na.rm=TRUE))

nunique.fun <- function(beta) length(unique(beta))

bias_n_fun <- function(beta_sim,beta_true,n=20){
  top_n <- order(abs(beta_sim),decreasing=TRUE)[1:n]
  beta_est_n <- beta_sim[top_n]
  beta_true_n <- beta_true[top_n] 
  mean(abs((beta_true_n-beta_est_n)/beta_true_n))
}

smc <- function(x,y){
  # x and y must be integer vectors with entries from 1 to number of groups
  # where number of groups may be different for each clustering
  # x= "true" reference grouping
  
  # pairs.x is a matrix where pairs.x[i,j]=1 if i and j are in the same cluster in x
  pairs.x <- Matrix(outer(X=x,Y=x,FUN=function(x,y) x==y))

  # pairs.y is a matrix where pairs.y[i,j]=1 if i and j are in the same cluster in y
  pairs.y <- Matrix(outer(X=y,Y=y,FUN=function(x,y) x==y))
  
  # when do x and y agree on which outcomes are paired together?
  pairs.xy <- pairs.x == pairs.y
  # remove the lower triangle (including diagonal)
  # pairs.xy[lower.tri(pairs.xy,diag=T)] <- FALSE
  diag(pairs.xy) <- FALSE
  pairs.count <- Matrix(TRUE,nrow=length(x),ncol=length(x))
  diag(pairs.count) <- FALSE
  
  # get average agreement for each group in x 
  num <- numeric(length=max(x))
  denom <- numeric(length=max(x))
  n <- length(x)
  for(i in 1:max(x)){
    n.i <- sum(x==i)
    n.ni <- n-n.i
    num[i] <- sum(pairs.xy[x==i,])
    denom[i] <- sum(pairs.count[x==i,])
  }
  smc.xy <- num/denom
  # # two numbers below should be the same
  # require(fossil)
  # rand.index(x,y)
  # sum(num)/(2*choose(n,2))
  return(smc.xy)
}

# purity <- function(x,y){
#   # x and y must be integer vectors with entries from 1 to number of groups
#   # where number of groups may be different for each clustering
#   # x is "true" set of clusters to compare to new clustering, y
#   p <- numeric(max(x))
#   which.y <- numeric(max(x))
#   for(i in 1:max(x)){
#     p.i <- sapply(1:max(y),function(z,x,y) mean(y[x==i]==z),
#                   x=x,y=y)
#     which.y.i <- (1:max(y))[which.max(p.i)]
#     p[i] <- p.i[which.y.i]
#     which.y[i] <- which.y.i
#   }
#   names(p) <- as.character(1:max(x))
#   df <- data.frame(x=1:max(x),y=which.y,p=p)
#   return(df)
# }

# jidx <- function(x,y){
#   # x and y must be integer vectors with entries from 1 to number of groups
#   # where number of groups may be different for each clustering
#   # x= "true" reference grouping
#   pairs.x <- Matrix(outer(X=x,Y=x,FUN=function(x,y) x==y))
#   pairs.x[lower.tri(pairs.x,diag = T)] <- FALSE # note: we keep TRUE on diagonal due to single outcome groups
#   pairs.y <- Matrix(outer(X=y,Y=y,FUN=function(x,y) x==y))
#   pairs.y[lower.tri(pairs.y,diag=T)] <- FALSE
#   pairs.num <- pairs.x & pairs.y # which pairs of outcomes are assigned to same group by both x and y
#   pairs.denom <- pairs.x | pairs.y # which pairs of outcomes are assigned to same group by at least one of x and y
#   # ## Check: the two values below should be the same
#   # sum(pairs.num[lower.tri(pairs.num)])/(sum(pairs.denom[lower.tri(pairs.denom)]))
#   # cluster_similarity(x,y) # Jaccard index according to clusteval package
#   # ##
#   smc.xy <- numeric(length=max(x))
#   for(i in 1:max(x)){
#     smc.xy[i] <- sum(pairs.num[x==i,])/sum(pairs.denom[x==i,])
#   }
#   return(smc.xy)
# }