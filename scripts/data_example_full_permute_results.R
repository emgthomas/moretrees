# --------------------------------------------------------------------------------- #
# ------------------- Results: permuting case/control status ---------------------- #
# --------------------------------------------------------------------------------- #

# direc <- "../moretrees/" # path of the moretrees repository
direc <- "/Users/emt380/Documents/PhD_Papers/Air_pollution/R_code/MORETreeS/moretrees/"
# direc <- "/nfs/home/E/ethomas/shared_space/ci3_nsaph/Emma/R_code/moretrees"
setwd(direc)

if(!dir.exists("./data_example_results")) dir.create("./data_example_results")

##### Load functions ######
source("scripts/VI_functions.R")
source("scripts/processing_functions.R")
require(igraph)

##### Read in permutations #####

nperm <- 10
nrestarts <- 3

for(perm in 1:nperm){
  for(restart in 1:nrestarts){
    file <- paste0("./data_example_results/data_example_full_perm",perm,"_restart",restart,".Rdata")
    if(file.exists(file)){
      load(file)
      print(unique(out_ss$moretrees_est))
    }
  }
}


