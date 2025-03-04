library(parallel) #for using multiple cores (mclapply for mac and linux, parLapply for windows)
# library(randomForest)
# library(lhs)

# load function
source("funcs/00_model_cluster.R")

# load in data
# parm_mat <- readRDS("dat/mat_var_0131.Rdata")
parm_mat <- readRDS("dat/mat_var_0211.Rdata")
head(parm_mat)
# will need 500 iterations of each parmeter combo
N <- dim(parm_mat)[1]
# print(N)
rep_num <- 500 # need to run 3 times, change to 333 for next run sets
rep_parm <- matrix(rep(t(parm_mat), rep_num), ncol = ncol(parm_mat), byrow = T)

# make a csv file to write to
result_dat <- paste0("dat/sim_results_", format(Sys.time(), "%m%d"), ".csv")
file.create(result_dat)
val_names <- c("extinct,pop_drop20, pop_drop50, pop_drop80, r_allele_peak15, r_allele_peak45, r_allele_peak75, final_r_allele, final_pop_size, final_inf_prev, max_r_allele, time_max_r_allele, max_inf_prev, time_last_zero_inf, min_pop, time_min_pop, at_K95, firstK95, r_ts_d0, parm_number \n")
cat(val_names, file = result_dat, append = TRUE)

# run functions across parm_mat

# run Monte Carlo simulations over parameterizations
print("starting simulations")

sim_results <- mclapply(1:(N*rep_num), function(i){
  # if (i %% 200000 == 0) {print(i)} # just to keep track of where we are?
  case_vect <- rep_parm[i, c(1:34, 36)] # skip density dependence only, since that manipulates parameters when making the matrix
  tmp <- cluster_run(case_vect) 
  cat(paste(shQuote(tmp, type="cmd"), collapse=","), file=result_dat, append=TRUE) # save the data
  cat("\n", file=result_dat, append=TRUE) # add a line
}, mc.cores = parallel::detectCores()) 

print("simulations complete")
# print(length(sim_results)) #sim_results should be a list of 1e6 lists.... sorry
# should all be in dat/sim_results.csv hopefully!

# convert to df
# sim_dat <- as.data.frame(matrix(unlist(sim_results), ncol = 20, byrow = T))
# 
# colnames(sim_dat) <- c("extinct",
#                        "pop_drop20", "pop_drop50", "pop_drop80",
#                        "r_allele_peak15", "r_allele_peak45", "r_allele_peak75",
#                        "final_r_allele", "final_pop_size", "final_inf_prev",
#                        "max_r_allele", "time_max_r_allele",
#                        "max_inf_prev", "time_last_zero_inf",
#                        "min_pop", "time_min_pop",
#                        "at_K95", "firstK95",
#                        "r_ts_d0", "parm_number")

# saveRDS(sim_dat, file = paste0("dat/sim_dat_", format(Sys.time(), "%m%d"), ".Rdata"))
# saveRDS(sim_results, file = paste0("dat/sim_dat_", format(Sys.time(), "%m%d"), ".Rdata"))
