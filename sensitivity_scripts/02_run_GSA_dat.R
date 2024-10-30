library(parallel) #for using multiple cores (mclapply for mac and linux, parLapply for windows)
# library(randomForest)

# load function
source("functions/model_cluster.R")
# load data
parm_mat <- readRDS("sensitivity/mat_var.rds")
# will need 500 iterations of each parmeter combo
N <- dim(parm_mat)[1]
rep_num <- 500
rep_parm <- matrix(rep(t(parm_mat), rep_num), ncol = ncol(parm_mat), byrow = T)

# run functions across parm_mat
# prep for using multiple clusters
cl <- detectCores()

# run Monte Carlo simulations over parameterizations
print("starting simulations")

sim_results <- mclapply(1:(N*rep_num), function(i){
  cluster_run(rep_parm[i, ]) 
}, mc.cores = n_cores) 

print("simulations complete")
print(length(sim_results)) #sim_results should be a list of 1e6 lists.... sorry

# convert to df
sim_dat <- as.data.frame(matrix(unlist(sim_results), ncol = 20, byrow = T))
colnames(sim_dat) <- c("extinct", 
                       "pop_drop20", "pop_drop50", "pop_drop80",
                       "r_allele_peak15", "r_allele_peak45", "r_allele_peak75",
                       "final_r_allele", "final_pop_size", "final_inf_prev",
                       "max_r_allele", "time_max_r_allele",
                       "max_inf_prev", "time_last_zero_inf",
                       "min_pop", "time_min_pop",
                       "at_K95", "firstK95",
                       "r_ts_d0", "parm_number")

# summarize by case
sens_dat_ext <- sim_dat %>% group_by(parm_number) %>% 
  summarise(`P(extinct)` = sum(extinct)/rep_num, 
            `P(pop_drop20)` = sum(pop_drop20)/rep_num,
            `P(pop_drop50)` = sum(pop_drop50)/rep_num,
            `P(pop_drop80)` = sum(pop_drop80)/rep_num, 
            avg_final_r = mean(final_r_allele),
            sd_final_r = sd(final_r_allele),
            avg_final_inf = mean(final_inf_prev),
            sd_final_inf = sd(final_inf_prev))
sens_dat_erE <- sim_dat %>% 
  filter(extinct == 0, pop_drop20 == 1, at_K95 == 1, r_allele_peak75 == 1) %>%
  group_by(parm_number) %>% 
  summarise(`P(ER_20_75)` = n()/rep_num)
tmp <- merge(sens_dat_ext, sens_dat_erE, by = "parm_number", all = T)
sens_dat_erM <- sim_dat %>% 
  filter(extinct == 0, pop_drop50 == 1, at_K95 == 1, r_allele_peak45 == 1) %>%
  group_by(parm_number) %>% 
  summarise(`P(ER_50_45)` = n()/rep_num)
tmp <- merge(tmp, sens_dat_erM, by = "parm_number", all = T)
sens_dat_erP <- sim_dat %>% 
  filter(extinct == 0, pop_drop80 == 1, at_K95 == 1, r_allele_peak15 == 1) %>%
  group_by(parm_number) %>% 
  summarise(`P(ER_80_15)` = n()/rep_num)
tmp <- merge(tmp, sens_dat_erP, by = "parm_number", all = T)
sens_dat_drE <- sim_dat %>% 
  filter(extinct == 0, pop_drop20 == 1, at_K95 == 1, r_allele_peak75 == 0) %>%
  group_by(parm_number) %>% 
  summarise(`P(DR_20_75)` = n()/rep_num)
tmp <- merge(tmp, sens_dat_drE, by = "parm_number", all = T)
sens_dat_drM <- sim_dat %>% 
  filter(extinct == 0, pop_drop50 == 1, at_K95 == 1, r_allele_peak45 == 0) %>%
  group_by(parm_number) %>% 
  summarise(`P(DR_50_45)` = n()/rep_num)
tmp <- merge(tmp, sens_dat_drM, by = "parm_number", all = T)
sens_dat_drP <- sim_dat %>% 
  filter(extinct == 0, pop_drop80 == 1, at_K95 == 1, r_allele_peak15 == 0) %>%
  group_by(parm_number) %>% 
  summarise(`P(DR_80_15)` = n()/rep_num)
tmp <- merge(tmp, sens_dat_drP, by = "parm_number", all = T)
sens_dat_maE <- sim_dat %>% 
  filter(extinct == 0, pop_drop80 == 0, at_K95 == 1, r_allele_peak75 == 0) %>%
  group_by(parm_number) %>% 
  summarise(`P(MA_75)` = n()/rep_num)
tmp <- merge(tmp, sens_dat_maE, by = "parm_number", all = T)
sens_dat_maM <- sim_dat %>% 
  filter(extinct == 0, pop_drop80 == 0, at_K95 == 1, r_allele_peak45 == 0) %>%
  group_by(parm_number) %>% 
  summarise(`P(MA_45)` = n()/rep_num)
tmp <- merge(tmp, sens_dat_maM, by = "parm_number", all = T)
sens_dat_maP <- sim_dat %>% 
  filter(extinct == 0, pop_drop80 == 0, at_K95 == 1, r_allele_peak15 == 0) %>%
  group_by(parm_number) %>% 
  summarise(`P(MA_15)` = n()/rep_num)
tmp <- merge(tmp, sens_dat_maP, by = "parm_number", all = T)
# note NA's are really 0s
GSA_summary_dat <- tmp

save(GSA_summary_dat, file = "sensitivity/GSA_summary_dat.Rdata")
