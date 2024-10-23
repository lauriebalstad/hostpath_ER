library(parallel)
load("sensitivity/mat_var.rds")
source("functions/model_cluster.R")

### set up data: need to rep 500 times each parameter combination
parm_mat <- matrix(rep(t(mat_var), 500), ncol = ncol(mat_var), byrow = TRUE)
m <- dim(parm_mat)[1] 

### run simulations -- will have 500*2000 = 1e6 simulations
sim_results <- mclapply(1:m, 
                        function(i){cluster_run(parm_mat[i, 1:29], # parm_vect
                                                floor(parm_mat[i, 30]/parm_mat[i, 29]), # ngens
                                                parm_mat[i, 31])}, # parm_id
                        mc.cores = parallel::detectCores()-1) 

### deal with outputs (list of lists)