# ------------------------------------------------------- #
# ------- Create folds for cross-validation ------------- #
# ------------------------------------------------------- #

# direc <- "/nfs/home/E/ethomas/shared_space/ci3_nsaph/Emma/R_code/MORETreeS/"
# direc <- "/Users/emt380/Documents/PhD_Papers/Air_pollution/R_code/MORETreeS/"
direc <- "../moretrees/" # path of the moretrees repository
setwd(direc)

# Load data
#load(file="/nfs/home/E/ethomas/shared_space/ci3_nsaph/Emma/Data/Case_crossover_data/moretrees_CC_data.Rdata")
load(file="data/moretrees_CC_data.Rdata")

# Set seed
set.seed(943860)

# Split data into nfolds "folds" stratified by each outcome v
folds <- list()
for(v in 1:pL){
  n.v <- Y[v]
  folds[[v]] <- sample(1:nfolds,n.v,replace=T)
}

save(folds,file="./data/cv_folds.Rdata")