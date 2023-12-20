# this is getting all the parameters and such together

library(dplyr) # organizing data
library(randomForest) # sensitivity
library(lhs) # sensitivity
library(ggplot2)

# drawing parameters for sensitivity
# sample size
N <- 200 # testing....
# parameters -- fix ngen to 45, number of sims per situation to 500
# remember that parameters are in rates
params <- c("comp_str", "event_order", # structural
            "n_RR", "n_WR", "n_WW", "init_I", # set ups
            "b_RR", "b_WR", "b_WW", "b_sd", # transmission
            "m_SRR", "m_SWR", "m_SWW", "m_IRR", "m_IWR", "m_IWW", "m_Ssd", "m_Isd", # mortality
            "r_RR", "r_WR", "r_WW", "r_sd", # recovery
            "l_RR", "l_WR", "l_WW", "l_sd", # reproduction
            "K", "K_sd") # carrying capacity

mat <- randomLHS(n = N, k = length(params))# need to update the matrices so that the right sorts of things are in the matrix
mat[,1] <- ceiling(qunif(mat[,1], 0, 3))
mat[,2] <- ceiling(qunif(mat[,2], 0, 6))
mat[,3] <- ceiling(qunif(mat[,3], 0, 10)) # has to be positive whole numbers
mat[,4] <- ceiling(qunif(mat[,4], 0, 10)) # has to be positive whole numbers
mat[,5] <- ceiling(qunif(mat[,5], 0, 10)) # has to be positive whole numbers
# let init_I be any number 0-1 --> okay
mat[,9] <- qnorm(mat[,9], 1, 0.2) # transform into rate with mean 1 and sd 0.1
mat[,7] <- mat[,9]*mat[,7] # b_RR as percent of b_WW value
mat[,8] <- qunif(mat[,8],min=mat[,7],max=mat[,9]) # set b_WR between
# leave sd
mat[,13] <- qnorm(mat[,13], 0.1, 0.02) # transform into rate with mean 0.1 and sd 0.02
mat[,11] <- mat[,13]*mat[,11] # b_RR as percent of m_WW value
mat[,12] <- qunif(mat[,12],min=mat[,11],max=mat[,13]) # set m_WR between
# m_IXX is additive
mat[,16] <- qnorm(mat[,16], 0.9, 0.05) # transform into rate with mean 0.9 and sd 0.05
mat[,14] <- mat[,16]*mat[,14] # b_RR as percent of m_WW value
mat[,15] <- qunif(mat[,15],min=mat[,14],max=mat[,16]) # set m_WR between
# leave both sd
mat[,21] <- qnorm(mat[,21], 1, 0.1) # transform into rate with mean 1 and sd 0.1
mat[,19] <- mat[,21]*mat[,19] # r_RR as percent of b_WW value
mat[,20] <- qunif(mat[,20],min=mat[,19],max=mat[,21]) # set r_WR between
# leave sd
mat[,25] <- qnorm(mat[,25], 4, 0.1) # highest fitness for WWs
mat[,23] <- mat[,25]*mat[,23] # l_RR as percent of l_WW value
mat[,24] <- qunif(mat[,24],min=mat[,23],max=mat[,25]) # set l_WR between
# leave sd
mat[,27] <- qnorm(mat[,27], 180, 20) # keep pretty flat
mat[,28] <- qnorm(mat[,28], 10, 1) # keep wide-ish
# check for zeros
mat[which(mat<0)] = 0

# run the matrix
sens_dat_1 <- mapply(run_gens, mat[1:100,1], mat[1:100,2],
                     mat[1:100,3], mat[1:100,4], mat[1:100,5], mat[1:100,6], 
                     mat[1:100,7], mat[1:100,8], mat[1:100,9], mat[1:100,10], 
                     mat[1:100,11], mat[1:100,12], mat[1:100,13], mat[1:100,14], mat[1:100,15], mat[1:100,16], mat[1:100,17], mat[1:100,18],
                     mat[1:100,19], mat[1:100,20], mat[1:100,21], mat[1:100,22], 
                     mat[1:100,23], mat[1:100,24], mat[1:100,25], mat[1:100,26], 
                     mat[1:100,27], mat[1:100,28], 
                     45, 400)
sens_dat_2 <- mapply(run_gens, mat[101:200,1], mat[101:200,2],
                     mat[101:200,3], mat[101:200,4], mat[101:200,5], mat[101:200,6], 
                     mat[101:200,7], mat[101:200,8], mat[101:200,9], mat[101:200,10], 
                     mat[101:200,11], mat[101:200,12], mat[101:200,13], mat[101:200,14], mat[101:200,15], mat[101:200,16], mat[101:200,17], mat[101:200,18],
                     mat[101:200,19], mat[101:200,20], mat[101:200,21], mat[101:200,22], 
                     mat[101:200,23], mat[101:200,24], mat[101:200,25], mat[101:200,26], 
                     mat[101:200,27], mat[101:200,28], 
                     45, 400)

sens_dat <- rbind(sens_dat_1, sens_dat_2)

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






