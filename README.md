# moretrees
Code for implementing Multi-outcome regression with Tree-Structured Shrinkage (MOReTreeS), as described in the manuscript of the same title.

This repository contains the following directories. The most important directory is scripts, which contains R files for implementing the model.

## data
Should contain case-crossover study data examining the effect of short-term exposure to PM2.5 on hospitalizations for cardiovascular disease (CVD) among Medicare enrollees. For data privacy reasons we are not able to make these data publicly available.

## data_example_results
Stores results for the data application.

## figures_and_tables
Contains pdf and tex files respectively for the figures and tables shown in the manuscript.

## scripts
Code for fitting spike and slab MOReTreeS (ssMOReTreeS) to data via a variational inference algorithm. Also includes scripts for running the simulations and data application shown in the manuscript, as well as for reproducing relevant tables and figures.

## simulation_inputs
Contains some R objects required as inputs to the simulation study. This include the tree of outcomes stored as an igraph object.

## simulation_results
Stores all output from the simulations.

