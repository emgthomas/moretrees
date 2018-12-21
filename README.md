# moretrees
Code for implementing Multi-outcome regression with Tree-Structured Shrinkage (MOReTreeS)

This repository contains the following files.

## VI_functions.R
Functions necessary for fitting spike & slab MOReTreeS (ssMOReTreeS) model to pair matched case-control data, including:

### adhoc_collapsing(Z,Y,pL,groups)
Computes log odds ratios for pair-matched case control data via conditional maximum likelihood and various adhoc collasping strategies specified by groups.

### initial_node_coeffs(Z,Y,uncollapsed,p,pL,leaf.descendants,ancestors)
Computes appropriate initial values for the variational parameter mu_gamma based on conditional logistic regression estimates of the log odds ratio for different levels of collapsing according to the tree.

### g_fun(eta)
Evaluates the function g (defined on pg 5 of the Supplementary Material)

### g_fun.vec(eta)
Vectorized version of g_fun.

### ELBO.fun_ss(Y,Z,p,pL,n,ancestors,VI_params,hyperparams,ELBO_old,tol,update_hyper)
Computes the evidence lower bound (ELBO) for ssMOReTreeS.

### VI_step_ss(ELBO,VI_params,hyperparams,Z,Y,n,p,pL,ancestors,leaf.descendants,update_hyper,tol)
Performs updates for variational parameters and hyperparameters in the ssMOReTreeS model.

### VI_binary_ss(Z,Y,n,p,pL,ancestors,leaf.descendants,cutoff,mu_gamma_init,tol,m.max,m.print,more,update_hyper,update_hyper_freq=)
Runs variational inference algorithm until convergence.

## data_example_full.R
Code for fitting ssMOReTreeS to case-crossover study data examining the effect of short-term exposure to PM2.5 on hospitalizations for cardiovascular disease (CVD) among Medicare enrollees.

## data_example_cv.R
Code for 10 fold cross-validation comparing predictive performance of ssMOReTreeS to various adhoc collapsing strategies with maximum likelihood fit. Models are fit to the same dataset as descried in data_example_full.R.

