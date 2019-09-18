# moretrees
Code for implementing Multi-outcome regression with Tree-Structured Shrinkage (MOReTreeS), as described in the manuscript of the same title.

This repository contains the following directories. The most important directory is scripts, which contains R files for implementing the model.

Data analyses cannot be fully reproduced as Medicare data cannot be made publicly available. However, our simulation study (Section 4 of the manuscript) can be reproduced via the scripts/simulations.R script.

## data
Should contain case-crossover study data examining the effect of short-term exposure to PM2.5 on hospitalizations for cardiovascular disease (CVD) among Medicare enrollees. For data privacy reasons we are not able to make these data publicly available.

## data_example_results
Stores results for the data application.

## figures_and_tables
Contains pdf and tex files respectively for the figures and tables shown in the manuscript.

## output
Outputs produced by the Condor batch processing system used to submit larger jobs to a cluster.

## scripts
Code for fitting spike and slab MOReTreeS (ssMOReTreeS) to data via a variational inference algorithm. Also includes scripts for running the simulations and data application shown in the manuscript, as well as for reproducing relevant tables and figures.

## simulation_inputs
Contains some R objects required as inputs to the simulation study and data example. This include the tree of outcomes stored as an igraph object.

## simulation_results
Stores all output from the simulations.

## submit_files
Submit files used to submit R code in ../scripts directory as jobs to a cluster via the Condor batch processing system.
