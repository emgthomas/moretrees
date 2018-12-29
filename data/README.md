# data

This directory contains case-crossover study data examining the effect of short-term exposure to PM2.5 on hospitalizations for cardiovascular disease (CVD) among Medicare enrollees.

The following files should be in this directory. However, for data privacy reasons we are not able to make these data publicly available.

## moretrees_CC_data.Rdata
Contains two key objects:
Z: a list of length p_L (number of outcomes). Each element of the list contains a numeric vector of length n_v (number of case-crossover pairs for each outcome), where each element is the difference in PM2.5 (in micrograms per cubic meter) between the case day and the control day.
Y: a numeric vector of lenght p_L, containing the number of cases associated with each outcome (n_v).

## cv_folds.Rdata
For use in cross-validation to examine the predictive validity of ssMOReTreeS in the data example. This file can be created by the script data_example_cv_folds.R. Contains one object:
folds: a list of length p_L. Each element of the list contains a numeric vector of length n_v (number of case-crossover pairs for each outcome), where each element is a number from 1 to 10 indicating the fold that the observation belongs to.
