# --------------------------------------------------------------------------------- #
# --------------------- Data analysis: creating permutations ---------------------- #
# --------------------------------------------------------------------------------- #

# direc <- "../moretrees/" # path of the moretrees repository
# direc <- "/Users/emt380/Documents/PhD_Papers/Air_pollution/R_code/MORETreeS/moretrees/"
direc <- "/nfs/home/E/ethomas/shared_space/ci3_nsaph/Emma/R_code/moretrees"
setwd(direc)

### Load ICD9 tree
load("./simulation_inputs/inputs.Rdata")

### Create permutations of dataset to test robustness of results to tree structure
set.seed(5685484)
n.perm <- 10 # number of permutations
permutations <- matrix(nrow=pL,ncol=10)
for(i in 1:10){
  permutations[,i] <- sample(1:pL,pL,replace=F)
}

### Save permutations
save(permutations,file="./simulation_inputs/permutations.Rdata")