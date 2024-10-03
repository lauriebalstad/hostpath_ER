# this is getting all the parameters and such together

library(dplyr) # organizing data
library(randomForest) # sensitivity
library(ggplot2)
library(parallel)
library(foreach)
library(doParallel)

# drawing parameters for sensitivity
# sample size
N <- 2000 # testing....
# parameters -- fix ngen to 45, number of sims per situation to 500
# remember that parameters are in rates
params <- c("comp_str", "event_order", # structural
            "n_RR", "n_WR", "n_WW", "init_I", # set ups
            "b_RR", "b_WR", "b_WW", "b_sd", # transmission
            "m_SRR", "m_SWR", "m_SWW", "m_IRR", "m_IWR", "m_IWW", "m_Ssd", "m_Isd", # mortality
            "r_RR", "r_WR", "r_WW", "r_sd", # recovery
            "l_RR", "l_WR", "l_WW", "l_sd", # reproduction
            "K", "K_sd") # carrying capacity

# # make cluster for running parallel -- helps things run more quickly (run time < 10% non-parallel)
cl <- parallel::makeCluster(detectCores())
doParallel::registerDoParallel(cl)

# run the matrix a bunch of times
# break matrix into rows of 20
for (i in 1:200) {
  
  # name the matrix
  sub_mat_name <- paste0("mat_",i, sep = "")
  row_1 <- (i-1)*10 + 1
  row_fin <- i*10
  # print(c(sub_mat_name, row_1, row_fin)) 
  listed <- unlist(as.list(mat[row_1:row_fin,])) 
  sublisted <- split(listed, rep(1:10, 28)[1:280]) 
  # need to convert from matrix to list of rows 
  sub_mat <- mapply(run_gens, sublisted, # parm_vect, ordered
                    ngens = 40, bnum = 700) # 40 generations each simulation, 700 simulations per parameter combo
  # rows are outputs, columns are cases (parms combos)
  saveRDS(sub_mat, file = paste0("~/GitHub/evorescue_hostpath_str/sensitivity_data/",sub_mat_name,".rds", sep =""))
  
}
# also try running in parallel -- would need to use foreach loop
# https://cran.r-project.org/web/packages/foreach/vignettes/foreach.html
# normal one takes ~2 hours RIP. with within loop parallel, takes about 10 minutes :) all runs should take ~30 hours
# thinking that doing the 700 reps in parallel is the move? not sure if this or the i's is the more effective move

stopCluster(cl)

colnames(sens_dat) <- c(params)
input_mat <- mat
colnames(input_mat) <- c("p_extinct", "p_demo_recovery", "avg_time_K", "sd_time_K",
                        "avg_pop_size_tf", "sd_pop_size_tf", "avg_pop_size_given_NE", "sd_pop_size_given_NE",
                        "p_lost_I", "avg_time_lost_I", "sd_time_lost_I",
                        "avg_infect_class_I_tf", "sd_infect_class_I_tf", "avg_infect_class_I_given_NE", "sd_infect_class_I_given_NE",
                        "avg_freq_R_tf", "sd_freq_R_tf", "avg_max_R", "avg_time_max_R", "gens", "reps")

sens_forest <- randomForest(sens_dat$p_extinct~.,  
                            data = input_mat,  
                            importance = TRUE, 
                            proximity = TRUE) 

print(sens_forest)
plot(sens_forest)
VarImpPlot(sens_forest)




