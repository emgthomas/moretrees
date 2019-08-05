# --------------------------------------------------------------------------------- #
# --------------------- Data analysis: creating permutations ---------------------- #
# --------------------------------------------------------------------------------- #

# direc <- "../moretrees/" # path of the moretrees repository
direc <- "/Users/emt380/Documents/PhD_Papers/Air_pollution/R_code/MORETreeS/moretrees/"
# direc <- "/nfs/home/E/ethomas/shared_space/ci3_nsaph/Emma/R_code/moretrees"
setwd(direc)

### Create permutations of dataset to test robustness of results to tree structure
set.seed(4563426)
n.perm <- 10 # number of permutations
permutations <- sample(1:1E9,n.perm,replace=F)

### Save permutations
save(permutations,file="./simulation_inputs/permutations.Rdata")