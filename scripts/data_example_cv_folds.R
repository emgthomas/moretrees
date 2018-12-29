# ------------------------------------------------------- #
# ------- Create folds for cross-validation ------------- #
# ------------------------------------------------------- #

direc <- "../moretrees/" # path of the moretrees repository
setwd(direc)

# Load data
load(file="data/moretrees_CC_data.Rdata")

# Set seed
set.seed(943860)

# Split data into nfolds "folds" stratified by each outcome v
folds <- list()
for(v in 1:pL){
  n.v <- Y[v]
  folds[[v]] <- sample(1:nfolds,n.v,replace=T)
}

save(folds,file="data/cv_folds.Rdata")