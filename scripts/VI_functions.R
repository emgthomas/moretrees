# --------------------------------------------------------------------------------- #
# -------- Variational Inference for MOReTreeS with Spike & Slab Prior ------------ #
# ------------------------- Case-crossover Likelihood ----------------------------- #
# --------------------------------------------------------------------------------- #

adhoc_collapsing <- function(Z,Y,pL,groups){
  
  ## Inputs ##
  
  # Z = list of length pL. Z[[v]] contains a numeric object of length n_v,
  # where v is the outcome and n_v is the number of observations (cases) corresponding
  # to outcome v. The entries of Z[[v]] are the difference in exposure values between
  # the case day and control day.
  
  # Y = numeric object of length pL. Y[v] = n_v, the number of observations (cases) corresponding
  # to outcome v.
  
  # pL = number of distinct outcomes (leaves) in the outcome tree
  
  # groups = matrix with pL rows and g columns, where g is the number of different
  # adhoc collapsing strategies to examine. Each column indicates the group to which each
  # outcome (indicated by the row) belongs for the corresponding collapsing strategy.
  
  ## Outputs ##
  
  # coeffs = a numeric matrix of size pL x g containing the estimated log odds ratios for
  # each outcome (row) and adhoc collapsing strategy (column). Note that the estimates will
  # be identical for all outcomes allocated to he same group.
  
  ## Code ##
  
  warn <- getOption("warn") # Suppress warnings due to glm.fit
  options(warn=-1)
  
  # Compute individual effect estimates for each outcome (no collapsing)
  uncollapsed <- numeric(length=pL)
  for(v in 1:pL){
    if(Y[v]==0){
      uncollapsed[v] <- 0
    } else {
      uncollapsed[v] <- glm(rep(1,Y[v]) ~ 0 + Z[[v]],family="binomial")$coefficients[1]
    }
  }
  
  # Compute grouped estimates, with groupings specified by the adhoc collapsing strategies
  # contained in the matrix groups
  beta.groups <- data.frame(matrix(nrow=nrow(groups),ncol=ncol(groups)))
  names(beta.groups) <- names(groups)
  for(i in 1:ncol(groups)){ # Loop through the adhoc collapsing strategies
    for(gr in levels(as.factor(groups[,i]))){ # Loop through the groups for the current strategy
      which.dat <- groups[,i]==gr
      Y.n <- sum(Y[which.dat])
      if(Y.n == 0){ # If the sample size for this group is zero, set the effect size to zero
        beta.groups[which.dat,i] <- 0
      } else { # If the sample size is non-zero, estimate the effect using conditional logistic regression
        beta_ml <- glm(rep(1,Y.n) ~ 0 + unlist(Z[which.dat]),family="binomial")$coefficients[1]
        beta.groups[which.dat,i] <- rep(beta_ml,sum(which.dat))
      }
    }
  }
  coeffs <- data.frame(uncollapsed=uncollapsed,beta.groups)
  
  options(warn=warn)
  
  return(coeffs)
}

initial_node_coeffs <- function(Z,Y,uncollapsed,p,pL,leaf.descendants,ancestors){
  
  ## Inputs ##
  
  # Z = list of length pL. Z[[v]] contains a numeric object of length n_v,
  # where v is the outcome and n_v is the number of observations (cases) corresponding
  # to outcome v. The entries of Z[[v]] are the difference in exposure values between
  # the case day and control day.
  
  # Y = numeric object of length pL. Y[v] = n_v, the number of observations (cases) corresponding
  # to outcome v.
  
  # uncollapsed = numeric object of length pL. uncollapsed[v] is the log odds ratio for the effect of
  # the exposure on outcome v estimated via standard maximum likelhood for the conditional logistic
  # regression model.
  
  # p = total number of nodes in the outcome tree (including leaves and internal nodes)
  
  # pL = number of distinct outcomes (leaves) in the outcome tree
  
  # leaf.descendants = list of length p, where leaf.descendants[[v]] is a list of integers indicating
  # the descendants of node v that are leaves (including v itself, if v is a leaf)
  
  # ancestors = list of length p, where ancestors[[p]] is a list of integers indicating the ancestors
  # of v on the tree (this list always includes v itself)
  
  ## Outputs ##
  
  # mu_gamma_init = a numeric object of length p, where mu_gamma_init[p] is an initial value for the
  # ssMOReTreeS model mu_gamma parameter.
  
  ## Code ##
  
  warn <- getOption("warn") # Suppress warnings due to glm.fit
  options(warn=-1)
  
  # Estimate the effects for different levels of collasping
  beta_init <- numeric(length=p)
  beta_init[(p-pL+1):p] <- uncollapsed
  for(v in 1:(p-pL)){
    u.vec <- leaf.descendants[[v]]-(p-pL)
    Y.n <- sum(Y[u.vec])
    if(Y.n == 0){
      beta_init[v] <- 0
    } else {
      beta_init[v] <- glm(rep(1,Y.n) ~ 0 + unlist(Z[u.vec]),family="binomial")$coefficients[1]
    }
  }
  
  # Find the difference in the above estimates between each node and its parent
  # These will be used as initial values for mu_gamma
  mu_gamma_init <- numeric(length=p)
  mu_gamma_init[1] <- beta_init[1]
  for(v in 2:p){
    if(beta_init[v] == 0){
      mu_gamma_init[v] <- 0
    } else {
      anc <- ancestors[[v]]
      parent <- anc[2]
      mu_gamma_init[v] <- beta_init[v] - beta_init[parent]
    }
  }
  
  options(warn=warn)
  
  return(mu_gamma_init)
}

g_fun <- function(eta){
  
  ## Inputs ##
  
  # eta = a variational parameter
  
  ## Outputs ##
  
  # g(eta), a function of eta
  
  ## Code ##
  
  (1/(2*eta))*(1/(1+exp(-eta))-0.5)
}

g_fun.vec <- function(eta){
  
  ## Inputs ##
  
  # eta = a numeric vector
  
  ## Outputs ##
  
  # g(eta), a numeric vector containg the values of the function g evaluated at each element of eta
  
  ## Code ##
  
  if(length(eta)==0){
    return(numeric(0))
  } else {
    y <- rep(1/8,length(eta))
    which.nonzero <- eta != 0
    y[which.nonzero] <- g_fun(eta[which.nonzero])
    return(y)
  }
}

log1p.exp <- function(x){
  
  ## Inputs ##
  
  # x = a numeric object of length 1
  
  ## Outputs ##
  
  # Returns the value log(1 + e^x)
  
  ## Code ##
  
  if(x > 20){
    return(x)
  } else {
    return(log1p(exp(x)))
  }
}

log1p.exp.vec <- function(x){
  
  ## Inputs ##
  
  # x = a numeric object of length > 1
  
  ## Outputs ##
  
  # Returns the value log(1 + e^x) for each element of x
  
  ## Code ##
  
  y <- x
  which.small <- x <= 20
  y[which.small] <- log1p(exp(x[which.small]))
  return(y)
}

loglogit <- function(x){
  
  ## Inputs ##
  
  # x = a numeric object of length > 1
  
  ## Outputs ##
  
  # Returns the value -1*log(1 + 1/e^(-x)) for each element of x
  
  ## Code ##
  
  -log1p.exp.vec(-x)
}

ELBO.fun_ss <- function(Y,Z,p,pL,n,ancestors,VI_params,hyperparams,ELBO_old=1,tol=1E-16,update_hyper=F){
  
  ## Inputs ##
  
  # Y = numeric object of length pL. Y[v] = n_v, the number of observations (cases) corresponding
  # to outcome v.'
  
  # Z = list of length pL. Z[[v]] contains a numeric object of length n_v,
  # where v is the outcome and n_v is the number of observations (cases) corresponding
  # to outcome v. The entries of Z[[v]] are the difference in exposure values between
  # the case day and control day.
  
  # p = total number of nodes in the outcome tree (including leaves and internal nodes)
  
  # pL = number of distinct outcomes (leaves) in the outcome tree
  
  # n = sum(Y), the total number of cases
  
  # VI_params = a list containing the variational parameters. Specifically, the list contains 
  # the named elements mu_gamma, sigma2_gamma and u_s, all numeric vectors of length p.
  
  # hyperparams = a list containing the hyperparameters. Specifically, the list contains the named
  # elements tau (numeric length 1), rho (numeric length 1) and eta (list of length p, where eta[[v]] is
  # a numeric of length Y[v]). eta is not actually a hyperparameter but is updated at the same time.
  
  # ancestors = list of length p, where ancestors[[p]] is a list of integers indicating the ancestors
  # of v on the tree (this list always includes v itself)
  
  # ELBO_old = the previous value of the evidence lower bound (ELBO)
  
  # tol = tolerance for converegence of the VI algorithm (difference between newly computed ELBO
  # and ELBO_old required to stop the algorithm)
  
  # update_hyper = boolean object. If update_hyper = TRUE, hyperparameters (rho and tau) will be updated
  # and returned as part of this function. If FALSE, only eta will be updated.
  
  ## Outputs ##
  
  # A list with the following named elements:
  
  # ELBO = value of the ELBO based on current variatonal parameters and hyperparameters
  
  # hyperparams = updated list of hyperparameters
  
  ## Code ##

  # Calculate some things needed below
  pi_s <- 1/(1+exp(-VI_params$u_s))
  sgamma_ <- VI_params$mu_gamma*pi_s
  gamma2_ <- (VI_params$mu_gamma^2 + VI_params$sigma2_gamma)*pi_s + hyperparams$tau*(1-pi_s)
  sgamma2_ <- (VI_params$mu_gamma^2 + VI_params$sigma2_gamma)*pi_s
  sgamma_mat <- sgamma_ %*% t(sgamma_)
  diag(sgamma_mat) <- 0
  eta.unlist <- unlist(hyperparams$eta)
  g.eta <- g_fun.vec(eta.unlist)
  
  # Compute ELBO line by line
  line1_vec <- numeric(pL)
  eta2_update_list <- sapply(Y,numeric)
  for(v in 1:pL){
    u.vec <- ancestors[[p-pL+v]]
    line1_vec[v] <- 0.5*sum(sgamma_[u.vec])*sum(Z[[v]])
    eta2_update_list[[v]] <- Z[[v]]^2*(sum(sgamma2_[u.vec]) + sum(sgamma_mat[u.vec,u.vec]))
  }
  
  line1 <- sum(line1_vec)
  
  line2 <- sum(loglogit(eta.unlist) - eta.unlist/2 + g.eta*eta.unlist^2)
  
  line3 <- -1*sum(g.eta*unlist(eta2_update_list))
  
  line4 <- log(hyperparams$rho)*sum(pi_s) + log(1-hyperparams$rho)*sum(1-pi_s) + sum(pi_s*0.5*log(2*pi*VI_params$sigma2_gamma))-sum(pi_s[pi_s!=0]*log(pi_s[pi_s!=0]))
  
  line5 <- 0.5*sum(pi_s)*(1- log(2*pi*hyperparams$tau)) - sum((1-pi_s[pi_s!=1])*log(1-pi_s[pi_s!=1])) - 0.5*sum(sgamma2_)/hyperparams$tau
  
  # Add them up!
  ELBO <- line1 + line2 + line3 + line4 + line5
  
  # Update eta
  hyperparams$eta <- sapply(eta2_update_list,sqrt)
  
  # Update hyperparameters
  if(update_hyper){
    # Only do this if have not yet reached convergence
    if(abs(ELBO-ELBO_old) > tol){
      hyperparams$tau <- mean(gamma2_)
      hyperparams$rho <- mean(pi_s)
    }
  }
  
  # Return
  return(list(ELBO=ELBO,hyperparams=hyperparams))
}

VI_step_ss <- function(ELBO,VI_params,hyperparams,Z,Y,n,p,pL,ancestors,leaf.descendants,update_hyper,tol=1E-16){
  
  ## Inputs ##
  
  # ELBO = most recently computed value of the evidence lower bound (ELBO)
  
  # VI_params = a list containing the variational parameters. Specifically, the list contains 
  # the named elements mu_gamma, sigma2_gamma and u_s, all numeric vectors of length p.
  
  # hyperparams = a list containing the hyperparameters. Specifically, the list contains the named
  # elements tau (numeric length 1), rho (numeric length 1) and eta (list of length p, where eta[[v]] is
  # a numeric of length Y[v]). eta is not actually a hyperparameter but is updated at the same time.
  
  # Z = list of length pL. Z[[v]] contains a numeric object of length n_v,
  # where v is the outcome and n_v is the number of observations (cases) corresponding
  # to outcome v. The entries of Z[[v]] are the difference in exposure values between
  # the case day and control day.
  
  # Y = numeric object of length pL. Y[v] = n_v, the number of observations (cases) corresponding
  # to outcome v.'
  
  # p = total number of nodes in the outcome tree (including leaves and internal nodes)
  
  # pL = number of distinct outcomes (leaves) in the outcome tree
  
  # n = sum(Y), the total number of cases
  
  # ancestors = list of length p, where ancestors[[p]] is a list of integers indicating the ancestors
  # of v on the tree (this list always includes v itself)
  
  # leaf.descendants = list of length p, where leaf.descendants[[v]] is a list of integers indicating
  # the descendants of node v that are leaves (including v itself, if v is a leaf)
  
  # update_hyper = boolean object. If update_hyper = TRUE, hyperparameters (rho and tau) will be updated
  # and returned as part of this function.
  
  # tol = tolerance for converegence of the VI algorithm (difference between newly computed ELBO
  # and ELBO_old required to stop the algorithm)
  
  ## Outputs ##
  
  # A list with the following named elements:
  
  # ELBO = updated value of the ELBO
  
  # epsilon = absolute difference between updated ELBO and previous ELBO
  
  # VI_params = updated list of variational parameters
  
  # hyperparams = updated list of hyperparameters
  
  ## Code ##
  
  ELBO_old <- ELBO
  
  # Compute current E_q[s[v]*gamma[v]] for each v
  sgamma_ <- 1/(1+exp(-VI_params$u_s))*VI_params$mu_gamma
  
  # ~.~ Updating spike & slab parameters ~.~ ##
  g.eta <- sapply(hyperparams$eta,g_fun.vec)
  for(v in 1:p){
    
    # Calculate A and B
    u.vec <- leaf.descendants[[v]]
    A <- 1/hyperparams$tau + 2*sum(unlist(Z[u.vec-(p-pL)])^2*unlist(g.eta[u.vec-(p-pL)]))
    B <- 0
    for(u in u.vec){
      w.vec <- ancestors[[u]]
      w.vec <- w.vec[w.vec != v]
      B <- B - 2*sum(sgamma_[w.vec])*sum(g.eta[[u-(p-pL)]]*Z[[u-(p-pL)]]^2) + 0.5*sum(unlist(Z[[u-(p-pL)]]))
    }
    
    # Update variational parameters for distribution of (s_v,gamma_v)
    VI_params$mu_gamma[v] <- B/A
    VI_params$sigma2_gamma[v] <- 1/A
    VI_params$u_s[v] <- -0.5*(log(A) + log(hyperparams$tau)) + log(hyperparams$rho) - log(1-hyperparams$rho) + (B^2)/(2*A)
    
    # Update sgamma_
    sgamma_[v] <- 1/(1+exp(-VI_params$u_s[v]))*VI_params$mu_gamma[v]
  }
  
  ## ~.~ Computing ELBO ~.~ ##
  
  ELBO <- ELBO.fun_ss(Y=Y,Z=Z,p=p,pL=pL,n=n,ancestors=ancestors,VI_params=VI_params,hyperparams=hyperparams,ELBO_old=ELBO_old,tol=tol,update_hyper=update_hyper)
  hyperparams <- ELBO$hyperparams
  ELBO <- ELBO$ELBO
  
  epsilon <- abs(ELBO - ELBO_old)
  
  # Return

  return(list(ELBO=ELBO,epsilon=epsilon,VI_params=VI_params,hyperparams=hyperparams))
  
}

VI_binary_ss <- function(Z,Y,n,p,pL,ancestors,leaf.descendants,cutoff=0.5,
                         mu_gamma_init=NULL,u_s_init=NULL,sigma2_gamma_init=NULL,tau_init=NULL,rho_init=NULL,
                         tol=1E-16,m.max=10000,m.print=m.max+1,more=FALSE,update_hyper=T,update_hyper_freq=10){
  
  ## Inputs ##
  
  # Z = list of length pL. Z[[v]] contains a numeric object of length n_v,
  # where v is the outcome and n_v is the number of observations (cases) corresponding
  # to outcome v. The entries of Z[[v]] are the difference in exposure values between
  # the case day and control day.
  
  # Y = numeric object of length pL. Y[v] = n_v, the number of observations (cases) corresponding
  # to outcome v.
  
  # n = sum(Y), the total number of cases
  
  # p = total number of nodes in the outcome tree (including leaves and internal nodes)
  
  # pL = number of distinct outcomes (leaves) in the outcome tree
  
  # ancestors = list of length p, where ancestors[[p]] is a list of integers indicating the ancestors
  # of v on the tree (this list always includes v itself)
  
  # leaf.descendants = list of length p, where leaf.descendants[[v]] is a list of integers indicating
  # the descendants of node v that are leaves (including v itself, if v is a leaf)
  
  # cutoff = minimum probability required for setting node selection parameters (s) to one when computing
  # approximate posterior effect estimates.
  
  # mu_gamma_init = numeric vector of length p containing initial values of the variational parameter mu_gamma
  
  # tol = tolerance in the evidence lower bound (ELBO) required for convergence of the VI algorithm
  
  # m.max = maximum number of time steps to run the algorithm.
  
  # m.print = the current time step will be printed out every m.print steps (default is no printing)
  
  # more = boolean object. If TRUE, the values of the variational parameters and hyperparameters at every
  # time step are returned. If FALSE (the default), only the final values are returned.
  
  # update_hyper = boolean object. If update_hyper = TRUE, hyperparameters (rho and tau) will be updated
  # as part of the algorithm. If FALSE, no hyperparameter updating will occur.
  
  # update_hyper_freq = integer dictating how often the hyperparameters will be updated. For example, if
  # update_hyper_freq = 10 (the default), the hyperparameters will be updated every 10 time steps. The default is 1.
  
  ## Outputs ##
  
  # If more = FALSE, a list with the following named elements. If more = TRUE, the corresponding objects are returned
  # as lists containing the values of the parameters at every time-step.
  
  # ELBO = final value of the ELBO
  
  # moretrees_est = numeric vector of length pL containing the collapsed effect estimates for each outcome
  
  # VI_params = final list of variational parameters
  
  # hyperparams = final list of hyperparameters
  
  # reached.max = Boolean. TRUE if the maximum number of time-steps was reached; else FALSE.
  
  ## Code ##
  
  # If initial values for mu_gamma were not supplied, initialize randomly
  if(is.null(mu_gamma_init)){
    mu_gamma_init <- rnorm(p,initialize_params$m_mu_init,initialize_params$sd_mu_init)
  }
  
  # If initial values for u_s_init were not supplied, initialize randomly
  if(is.null(u_s_init)){
    pi_init <- runif(p)
    u_s_init <- log(pi_init/(1-pi_init))
  }
  
  # If initial values for sigma2_gamma_init were not supplied, intialize randomly
  if(is.null(sigma2_gamma_init)){
    sigma2_gamma_init <- rgamma(p,1,1)
  }
  
  # Hyperparams 
  
  # If initial value for tau_init was not supplied, initialize randomly
  if(is.null(tau_init)){
    tau_init <- var(mu_gamma_init)
  }
  
  # If initial value for tau_init was not supplied, initialize randomly
  if(is.null(rho_init)){
    rho_init <- 0.5
  }
  
  # Randomly initialize eta (will be updated immediately)
  eta_init <- sapply(Y,rnorm)
  
  # Put initial values in list
  VI_params <- list(mu_gamma=mu_gamma_init,sigma2_gamma=sigma2_gamma_init,u_s=u_s_init)
  hyperparams <- list(eta=eta_init,tau=tau_init,rho=rho_init)
  
  # Get initial ELBO and first eta update
  ELBO <- ELBO.fun_ss(Y=Y,Z=Z,p=p,pL=pL,n=n,ancestors=ancestors,VI_params=VI_params,hyperparams=hyperparams,update_hyper=FALSE)
  hyperparams <- ELBO$hyperparams
  ELBO <- ELBO$ELBO
  
  # Initialize
  ELBO.M <- numeric(m.max)
  ELBO.M[1] <- ELBO
  hyperparams.M <- list(hyperparams)
  VI_params.M <- list(VI_params)
  epsilon <- 10
  m <- 1
  
  # Iterate
  while(epsilon > tol & m <= m.max){
    m <- m+1
    update_hyper_m <- (m %% update_hyper_freq == 0) & update_hyper
    if(more){
      VI_out <- VI_step_ss(ELBO=ELBO.M[m-1],VI_params=VI_params.M[[m-1]],hyperparams=hyperparams.M[[m-1]],Z=Z,Y=Y,n=n,p=p,pL=pL,
                           ancestors=ancestors,leaf.descendants=leaf.descendants,tol=tol,update_hyper=update_hyper_m)
      hyperparams.M[[m]] <- VI_out$hyperparams
      VI_params.M[[m]] <- VI_out$VI_params
    } else {
      VI_out <- VI_step_ss(ELBO=ELBO.M[m-1],VI_params=VI_params.M[[1]],hyperparams=hyperparams.M[[1]],Z=Z,Y=Y,n=n,p=p,pL=pL,
                           ancestors=ancestors,leaf.descendants=leaf.descendants,tol=tol,update_hyper=update_hyper_m)
      hyperparams.M <- list(VI_out$hyperparams)
      VI_params.M <- list(VI_out$VI_params)
    }
    ELBO.M[m] <- VI_out$ELBO
    epsilon <- VI_out$epsilon
    if(m %% m.print == 0) print(m)
  }
  print(paste("Converged after",m,"steps.",sep=" "))
  
  # Check if maximum number of time steps was reached
  reached.max <- (m >= m.max)
  
  ELBO.M <- ELBO.M[1:m]
  if(more){
    VI_params <- VI_params.M[[m]]
    hyperparams <- hyperparams.M[[m]]
  } else {
    VI_params <- VI_params.M[[1]]
    hyperparams <- hyperparams.M[[1]]
  }

  # Get final effect estimates for each outcome
  pi_s_final <- 1/(1+ exp(-VI_params$u_s))
  mu_gamma_final <- VI_params$mu_gamma
  beta_est <- numeric(length=pL)
  for(v in 1:pL){
    anc <- ancestors[[p-pL+v]]
    pi_v <- pi_s_final[anc]
    pi_v[pi_v<cutoff] <- 0
    pi_v[pi_v>=cutoff] <- 1
    beta_est[v] <- sum(VI_params$mu_gamma[anc]*pi_v)
  }
  
  # Return results
  
  if(more){
    return(list(ELBO=ELBO.M,moretrees_est=beta_est,VI_params=VI_params.M,hyperparams=hyperparams.M,reached.max=reached.max))
  }
  
  return(list(ELBO=ELBO.M[m],moretrees_est=beta_est,VI_params=as.data.frame(VI_params),hyperparams=c(rho=hyperparams$rho,tau=hyperparams$tau),reached.max=reached.max))
  
}
