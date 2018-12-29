# scripts

This directory contains R scripts for implementing the spike and slab multi-outcome regression with tree-structured shrinkage (ssMOReTreeS) model. Scripts are also provided for reproducing the simulation and data application results shown in the manuscript, including relevant tables and figures.

The following files are in this directory:

## data_example_full.R
Code for fitting ssMOReTreeS to case-crossover study data examining the effect of short-term exposure to PM2.5 on hospitalizations for cardiovascular disease (CVD) among Medicare enrollees. Results will be saved to a directory named data_example_results. Package dependencies: igraph, doParallel

## data_example_cv.R
Code for 10 fold cross-validation comparing predictive performance of ssMOReTreeS to various adhoc collapsing strategies with maximum likelihood fit. Models are fit to the same dataset as descried in data_example_full.R. Results will be saved to a directory named data_example_results.

## data_example_figures_and_tables.R
Code for reproducing tables and figures related to the data example. Figures and tables created will be saved to a directory named figures_and_tables. This script produces Figure 2, Figure 4, Figure S2, Figure S4, Table 1 and Supplementary Material Section C. Package dependencies: igraph, data.table, ggplot2, RColorBrewer, glue, collapsibleTree, circlize, icd, plotly, xtable,

## simulations.R
Code for running simulations described in the manuscript. Package dependencies: igraph, doParallel

## simulations_figures_and_tables.R
Code for reproducing tables and figures related to the simulations. Figures and tables created will be saved to a directory named figures_and_tables. This script produces Figure 3, Figure S1 Figure A2 and Table A1. Package dependencies: mclust, igraph, reshape2, ggplot2, xtable, glue, icd, collapsibleTree, circlize, plotly

## VI_functions.R
Functions necessary for fitting the ssMOReTreeS model to pair matched case-control data, including:

### adhoc_collapsing()
Computes log odds ratios for pair-matched case control data via conditional maximum likelihood and various adhoc collasping strategies specified by groups.

### initial_node_coeffs()
Computes appropriate initial values for the variational parameter mu_gamma based on conditional logistic regression estimates of the log odds ratio for different levels of collapsing according to the tree.

### g_fun()
Evaluates the function g (defined on pg 5 of the Supplementary Material)

### g_fun.vec()
Vectorized version of g_fun.

### log1p.exp(), log1p.exp.vec(), loglogit()
Functions for calculating the log of the logit of a value.

### ELBO.fun_ss()
Computes the evidence lower bound (ELBO) for ssMOReTreeS.

### VI_step_ss()
Performs updates for variational parameters and hyperparameters in the ssMOReTreeS model.

### VI_binary_ss()
Runs variational inference algorithm until convergence.

## processing_functions.R
Functions for processing the output of VI_binary_ss(). These are mainly for use in analyzing the simulation results and producing relevant tables and figures.

### indiv.beta.calc()
Function for computing individual level beta estimates, as opposed to collapsed estimates (See Section 3.3 of the manuscript)

### gamma.sim.fun()
Function to help with producing simulated creddible intervals for individual beta estimates

### indiv.beta.calc.ci()
Function for computing individual level beta estimates along with 95% credible intervals (See Section 3.3 of the manuscript)

### groups.calc.fun()
Function for computing extracting collapsed beta estimates along with 95% credible intervals (See Section 3.3 of the manuscript) and putting them into an igraph tree

### explainer.latex()
Function for expanding ICD9 codes to print nicely as a LaTeX table

### expand_groups_latex()
Function for printing out a group of ICD9 codes in LaTeX formatted according to their place in the ICD9 heirarchy.

### rmse.fun()
Computes root mean squared error for a vector estimator

### nunique.fun()
Computes number of unique elements in a vector

### bias_n_fun()
Computes the average fractional bias of the top n largest estimates in a list of estimates.


## collapsibleTreeNetwork_modified.R
Code for producing interactive trees. The main function in this file, collapsibleTreeNetwork2(), was adapted from the collapsibleTree package (version 0.1.6) created by Adeel Khan. The original version of this package is available via CRAN at the following link: https://CRAN.R-project.org/package=collapsibleTree. Package dependencies: icd

### explainer.zip(), expand_html(), expand_html_sims()
Functions for formatting the output of simulations and data example for display on interactive trees.

### collapsibleTreeNetwork2()
Produces interactive tree in html format.

