library(ggplot2)
library(dplyr)
library(tidyverse)

#-----READ IN DATA: 02/03-----
str_dat_03 <- readRDS("C:/Users/lauri/AppData/Roaming/MobaXterm/home/evorescue_hostpath_str/dat/str_fig_d2_0203.Rdata")
# convert to df
str_dtf_03 <- as.data.frame(matrix(unlist(str_dat_03), ncol = 20, byrow = T))
colnames(str_dtf_03) <- c("extinct", 
                       "pop_drop20", "pop_drop50", "pop_drop80",
                       "r_allele_peak15", "r_allele_peak45", "r_allele_peak75",
                       "final_r_allele", "final_pop_size", "final_inf_prev",
                       "max_r_allele", "time_max_r_allele",
                       "max_inf_prev", "time_last_zero_inf",
                       "min_pop", "time_min_pop",
                       "at_K95", "firstK95",
                       "r_ts_d0", "parm_number")

cb_dat_03 <- readRDS("C:/Users/lauri/AppData/Roaming/MobaXterm/home/evorescue_hostpath_str/dat/cb_fig_d2_0203.Rdata")
cb_dtf_03 <- as.data.frame(matrix(unlist(cb_dat_03), ncol = 20, byrow = T))
colnames(cb_dtf_03) <- c("extinct", 
                      "pop_drop20", "pop_drop50", "pop_drop80",
                      "r_allele_peak15", "r_allele_peak45", "r_allele_peak75",
                      "final_r_allele", "final_pop_size", "final_inf_prev",
                      "max_r_allele", "time_max_r_allele",
                      "max_inf_prev", "time_last_zero_inf",
                      "min_pop", "time_min_pop",
                      "at_K95", "firstK95",
                      "r_ts_d0", "parm_number")

dc_dat_03 <- readRDS("C:/Users/lauri/AppData/Roaming/MobaXterm/home/evorescue_hostpath_str/dat/dc_fig_d2_0203.Rdata")
dc_dtf_03 <- as.data.frame(matrix(unlist(dc_dat_03), ncol = 20, byrow = T))
colnames(dc_dtf_03) <- c("extinct", 
                      "pop_drop20", "pop_drop50", "pop_drop80",
                      "r_allele_peak15", "r_allele_peak45", "r_allele_peak75",
                      "final_r_allele", "final_pop_size", "final_inf_prev",
                      "max_r_allele", "time_max_r_allele",
                      "max_inf_prev", "time_last_zero_inf",
                      "min_pop", "time_min_pop",
                      "at_K95", "firstK95",
                      "r_ts_d0", "parm_number")

pm_dat_03 <- readRDS("C:/Users/lauri/AppData/Roaming/MobaXterm/home/evorescue_hostpath_str/dat/pm_fig_d2_0203.Rdata")
pm_dtf_03 <- as.data.frame(matrix(unlist(pm_dat_03), ncol = 20, byrow = T))
colnames(pm_dtf_03) <- c("extinct", 
                      "pop_drop20", "pop_drop50", "pop_drop80",
                      "r_allele_peak15", "r_allele_peak45", "r_allele_peak75",
                      "final_r_allele", "final_pop_size", "final_inf_prev",
                      "max_r_allele", "time_max_r_allele",
                      "max_inf_prev", "time_last_zero_inf",
                      "min_pop", "time_min_pop",
                      "at_K95", "firstK95",
                      "r_ts_d0", "parm_number")

#-----READ IN DATA: 02/05-----
str_dat_05 <- readRDS("C:/Users/lauri/AppData/Roaming/MobaXterm/home/evorescue_hostpath_str/dat/str_fig_d2_0205.Rdata")
# convert to df
str_dtf_05 <- as.data.frame(matrix(unlist(str_dat_05), ncol = 20, byrow = T))
colnames(str_dtf_05) <- c("extinct", 
                       "pop_drop20", "pop_drop50", "pop_drop80",
                       "r_allele_peak15", "r_allele_peak45", "r_allele_peak75",
                       "final_r_allele", "final_pop_size", "final_inf_prev",
                       "max_r_allele", "time_max_r_allele",
                       "max_inf_prev", "time_last_zero_inf",
                       "min_pop", "time_min_pop",
                       "at_K95", "firstK95",
                       "r_ts_d0", "parm_number")

cb_dat_05 <- readRDS("C:/Users/lauri/AppData/Roaming/MobaXterm/home/evorescue_hostpath_str/dat/cb_fig_d2_0205.Rdata")
cb_dtf_05 <- as.data.frame(matrix(unlist(cb_dat_05), ncol = 20, byrow = T))
colnames(cb_dtf_05) <- c("extinct", 
                      "pop_drop20", "pop_drop50", "pop_drop80",
                      "r_allele_peak15", "r_allele_peak45", "r_allele_peak75",
                      "final_r_allele", "final_pop_size", "final_inf_prev",
                      "max_r_allele", "time_max_r_allele",
                      "max_inf_prev", "time_last_zero_inf",
                      "min_pop", "time_min_pop",
                      "at_K95", "firstK95",
                      "r_ts_d0", "parm_number")

dc_dat_05 <- readRDS("C:/Users/lauri/AppData/Roaming/MobaXterm/home/evorescue_hostpath_str/dat/dc_fig_d2_0205.Rdata")
dc_dtf_05 <- as.data.frame(matrix(unlist(dc_dat_05), ncol = 20, byrow = T))
colnames(dc_dtf_05) <- c("extinct", 
                      "pop_drop20", "pop_drop50", "pop_drop80",
                      "r_allele_peak15", "r_allele_peak45", "r_allele_peak75",
                      "final_r_allele", "final_pop_size", "final_inf_prev",
                      "max_r_allele", "time_max_r_allele",
                      "max_inf_prev", "time_last_zero_inf",
                      "min_pop", "time_min_pop",
                      "at_K95", "firstK95",
                      "r_ts_d0", "parm_number")

pm_dat_05 <- readRDS("C:/Users/lauri/AppData/Roaming/MobaXterm/home/evorescue_hostpath_str/dat/pm_fig_d2_0205.Rdata")
pm_dtf_05 <- as.data.frame(matrix(unlist(pm_dat_05), ncol = 20, byrow = T))
colnames(pm_dtf_05) <- c("extinct", 
                      "pop_drop20", "pop_drop50", "pop_drop80",
                      "r_allele_peak15", "r_allele_peak45", "r_allele_peak75",
                      "final_r_allele", "final_pop_size", "final_inf_prev",
                      "max_r_allele", "time_max_r_allele",
                      "max_inf_prev", "time_last_zero_inf",
                      "min_pop", "time_min_pop",
                      "at_K95", "firstK95",
                      "r_ts_d0", "parm_number")

#-----CALC PROBS & MERGE LIKE DATA: EX-----
str_ex_03 <- str_dtf_03 %>% filter(extinct == 1) %>% group_by(parm_number) %>% summarise(`P(extinct) - 03` = n()/1000)
str_ex_05 <- str_dtf_05 %>% filter(extinct == 1) %>% group_by(parm_number) %>% summarise(`P(extinct) - 05` = n()/1000)
str_ex <- merge(str_ex_03, str_ex_05, by = "parm_number", all = T)
str_ex$diff <- ifelse(is.na(str_ex$`P(extinct) - 03`), 0, str_ex$`P(extinct) - 03`) - ifelse(is.na(str_ex$`P(extinct) - 05`), 0, str_ex$`P(extinct) - 05`)

pm_ex_03 <- pm_dtf_03 %>% filter(extinct == 1) %>% group_by(parm_number) %>% summarise(`P(extinct) - 03` = n()/1000)
pm_ex_05 <- pm_dtf_05 %>% filter(extinct == 1) %>% group_by(parm_number) %>% summarise(`P(extinct) - 05` = n()/1000)
pm_ex <- merge(pm_ex_03, pm_ex_05, by = "parm_number", all = T)
pm_ex$diff <- ifelse(is.na(pm_ex$`P(extinct) - 03`), 0, pm_ex$`P(extinct) - 03`) - ifelse(is.na(pm_ex$`P(extinct) - 05`), 0, pm_ex$`P(extinct) - 05`)

dc_ex_03 <- dc_dtf_03 %>% filter(extinct == 1) %>% group_by(parm_number) %>% summarise(`P(extinct) - 03` = n()/1000)
dc_ex_05 <- dc_dtf_05 %>% filter(extinct == 1) %>% group_by(parm_number) %>% summarise(`P(extinct) - 05` = n()/1000)
dc_ex <- merge(dc_ex_03, dc_ex_05, by = "parm_number", all = T)
dc_ex$diff <- ifelse(is.na(dc_ex$`P(extinct) - 03`), 0, dc_ex$`P(extinct) - 03`) - ifelse(is.na(dc_ex$`P(extinct) - 05`), 0, dc_ex$`P(extinct) - 05`)

cb_ex_03 <- cb_dtf_03 %>% filter(extinct == 1) %>% group_by(parm_number) %>% summarise(`P(extinct) - 03` = n()/1000)
cb_ex_05 <- cb_dtf_05 %>% filter(extinct == 1) %>% group_by(parm_number) %>% summarise(`P(extinct) - 05` = n()/1000)
cb_ex <- merge(cb_ex_03, cb_ex_05, by = "parm_number", all = T)
cb_ex$diff <- ifelse(is.na(cb_ex$`P(extinct) - 03`), 0, cb_ex$`P(extinct) - 03`) - ifelse(is.na(cb_ex$`P(extinct) - 05`), 0, cb_ex$`P(extinct) - 05`)

fig_ex <- rbind(str_ex, pm_ex, dc_ex, cb_ex)

#-----CALC PROBS & MERGE LIKE DATA: ER-----
str_er_03 <- str_dtf_03 %>% filter(extinct == 0 & pop_drop50 == 1 &  at_K95 == 1 & r_allele_peak45 == 1) %>% group_by(parm_number) %>% summarise(`P(ER) - 03` = n()/1000)
str_er_05 <- str_dtf_05 %>% filter(extinct == 0 & pop_drop50 == 1 &  at_K95 == 1 & r_allele_peak45 == 1) %>% group_by(parm_number) %>% summarise(`P(ER) - 05` = n()/1000)
str_er <- merge(str_er_03, str_er_05, by = "parm_number", all = T)
str_er$diff <- ifelse(is.na(str_er$`P(ER) - 03`), 0, str_er$`P(ER) - 03`) - ifelse(is.na(str_er$`P(ER) - 05`), 0, str_er$`P(ER) - 05`)

pm_er_03 <- pm_dtf_03 %>% filter(extinct == 0 & pop_drop50 == 1 &  at_K95 == 1 & r_allele_peak45 == 1) %>% group_by(parm_number) %>% summarise(`P(ER) - 03` = n()/1000)
pm_er_05 <- pm_dtf_05 %>% filter(extinct == 0 & pop_drop50 == 1 &  at_K95 == 1 & r_allele_peak45 == 1) %>% group_by(parm_number) %>% summarise(`P(ER) - 05` = n()/1000)
pm_er <- merge(pm_er_03, pm_er_05, by = "parm_number", all = T)
pm_er$diff <- ifelse(is.na(pm_er$`P(ER) - 03`), 0, pm_er$`P(ER) - 03`) - ifelse(is.na(pm_er$`P(ER) - 05`), 0, pm_er$`P(ER) - 05`)

dc_er_03 <- dc_dtf_03 %>% filter(extinct == 0 & pop_drop50 == 1 &  at_K95 == 1 & r_allele_peak45 == 1) %>% group_by(parm_number) %>% summarise(`P(ER) - 03` = n()/1000)
dc_er_05 <- dc_dtf_05 %>% filter(extinct == 0 & pop_drop50 == 1 &  at_K95 == 1 & r_allele_peak45 == 1) %>% group_by(parm_number) %>% summarise(`P(ER) - 05` = n()/1000)
dc_er <- merge(dc_er_03, dc_er_05, by = "parm_number", all = T)
dc_er$diff <- ifelse(is.na(dc_er$`P(ER) - 03`), 0, dc_er$`P(ER) - 03`) - ifelse(is.na(dc_er$`P(ER) - 05`), 0, dc_er$`P(ER) - 05`)

cb_er_03 <- cb_dtf_03 %>% filter(extinct == 0 & pop_drop50 == 1 &  at_K95 == 1 & r_allele_peak45 == 1) %>% group_by(parm_number) %>% summarise(`P(ER) - 03` = n()/1000)
cb_er_05 <- cb_dtf_05 %>% filter(extinct == 0 & pop_drop50 == 1 &  at_K95 == 1 & r_allele_peak45 == 1) %>% group_by(parm_number) %>% summarise(`P(ER) - 05` = n()/1000)
cb_er <- merge(cb_er_03, cb_er_05, by = "parm_number", all = T)
cb_er$diff <- ifelse(is.na(cb_er$`P(ER) - 03`), 0, cb_er$`P(ER) - 03`) - ifelse(is.na(cb_er$`P(ER) - 05`), 0, cb_er$`P(ER) - 05`)

fig_er <- rbind(str_er, pm_er, dc_er, cb_er)

#-----PLOT HISTS-----
# add column and pivot longer
fig_ex$metric <- rep("P(extinct)")
fig_er$metric <- rep("P(ER)")
fig_tot <- rbind(fig_ex[, c(1, 4, 5)], fig_er[, c(1, 4, 5)])
fig_long <- pivot_longer(fig_tot, cols = 2) 

# box plots?
comp_sim_box <- ggplot(fig_long, aes(x = value, fill = metric)) + 
  geom_boxplot(col = "black") + 
  scale_fill_manual(values = c("#ac1457", "#f1c4a2")) + 
  facet_wrap(~metric, nrow = 1) + 
  labs(x = "difference across two test simulations") + 
  theme_bw() + theme(legend.position = "none")

# histograms?
comp_sim_hist <- ggplot(fig_long, aes(x = value, fill = metric)) + 
  geom_histogram(col = "black", bins = 30) + 
  scale_fill_manual(values = c("#ac1457", "#f1c4a2")) + 
  facet_wrap(~metric, nrow = 1) + 
  labs(x = "difference across two test simulations", y = NULL) + 
  theme_bw() + theme(legend.position = "none", 
                     axis.text.x = element_text(angle = 45, hjust = 1))

# save figure
png("figs/figure_plot/comp_sim_hist_0203-05.png",height=85,width=170,res=400,units='mm')
print(comp_sim_hist)
dev.off()

quantile(fig_ex$diff, c(0.025, 0.975))
quantile(fig_er$diff, c(0.025, 0.975))