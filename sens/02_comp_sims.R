library(ggplot2)
library(dplyr)
library(tidyverse)

#-----READ IN DATA: 02/03-----
str_dat_A <- readRDS("dat/str_fig_0909.Rdata")
# convert to df
str_dtf_A <- as.data.frame(matrix(unlist(str_dat_A), ncol = 29, byrow = T))
colnames(str_dtf_A) <- c("extinct", 
                         "pop_drop20", "pop_drop50", "pop_drop80",
                         "r_allele_peak15", "r_allele_peak45", "r_allele_peak75",
                         "final_r_allele", "final_pop_size", "final_inf_prev",
                         "max_r_allele", "time_max_r_allele",
                         "max_inf_prev", "time_last_zero_inf",
                         "min_pop", "time_min_pop",
                         "first_20", "first_50", "first_80",
                         "last_20", "last_50", "last_80",
                         "tot_20", "tot_50", "tot_80",
                         "at_K95", "firstK95",
                         "r_ts_d0", "parm_number")

cb_dat_A <- readRDS("dat/cb_fig_0909.Rdata")
cb_dtf_A <- as.data.frame(matrix(unlist(cb_dat_A), ncol = 29, byrow = T))
colnames(cb_dtf_A) <- c("extinct", 
                        "pop_drop20", "pop_drop50", "pop_drop80",
                        "r_allele_peak15", "r_allele_peak45", "r_allele_peak75",
                        "final_r_allele", "final_pop_size", "final_inf_prev",
                        "max_r_allele", "time_max_r_allele",
                        "max_inf_prev", "time_last_zero_inf",
                        "min_pop", "time_min_pop",
                        "first_20", "first_50", "first_80",
                        "last_20", "last_50", "last_80",
                        "tot_20", "tot_50", "tot_80",
                        "at_K95", "firstK95",
                        "r_ts_d0", "parm_number")

dc_dat_A <- readRDS("dat/dc_fig_0910.Rdata")
dc_dtf_A <- as.data.frame(matrix(unlist(dc_dat_A), ncol = 29, byrow = T))
colnames(dc_dtf_A) <- c("extinct", 
                        "pop_drop20", "pop_drop50", "pop_drop80",
                        "r_allele_peak15", "r_allele_peak45", "r_allele_peak75",
                        "final_r_allele", "final_pop_size", "final_inf_prev",
                        "max_r_allele", "time_max_r_allele",
                        "max_inf_prev", "time_last_zero_inf",
                        "min_pop", "time_min_pop",
                        "first_20", "first_50", "first_80",
                        "last_20", "last_50", "last_80",
                        "tot_20", "tot_50", "tot_80",
                        "at_K95", "firstK95",
                        "r_ts_d0", "parm_number")

pm_dat_A <- readRDS("dat/pm_fig_0909.Rdata")
pm_dtf_A <- as.data.frame(matrix(unlist(pm_dat_A), ncol = 29, byrow = T))
colnames(pm_dtf_A) <- c("extinct", 
                        "pop_drop20", "pop_drop50", "pop_drop80",
                        "r_allele_peak15", "r_allele_peak45", "r_allele_peak75",
                        "final_r_allele", "final_pop_size", "final_inf_prev",
                        "max_r_allele", "time_max_r_allele",
                        "max_inf_prev", "time_last_zero_inf",
                        "min_pop", "time_min_pop",
                        "first_20", "first_50", "first_80",
                        "last_20", "last_50", "last_80",
                        "tot_20", "tot_50", "tot_80",
                        "at_K95", "firstK95",
                        "r_ts_d0", "parm_number")

#-----READ IN DATA: 02/05-----
str_dat_B <- readRDS("dat/str_fig_0910.Rdata")
# convert to df
str_dtf_B <- as.data.frame(matrix(unlist(str_dat_B), ncol = 29, byrow = T))
colnames(str_dtf_B) <- c("extinct", 
                         "pop_drop20", "pop_drop50", "pop_drop80",
                         "r_allele_peak15", "r_allele_peak45", "r_allele_peak75",
                         "final_r_allele", "final_pop_size", "final_inf_prev",
                         "max_r_allele", "time_max_r_allele",
                         "max_inf_prev", "time_last_zero_inf",
                         "min_pop", "time_min_pop",
                         "first_20", "first_50", "first_80",
                         "last_20", "last_50", "last_80",
                         "tot_20", "tot_50", "tot_80",
                         "at_K95", "firstK95",
                         "r_ts_d0", "parm_number")

cb_dat_B <- readRDS("dat/cb_fig_0910.Rdata")
cb_dtf_B <- as.data.frame(matrix(unlist(cb_dat_B), ncol = 29, byrow = T))
colnames(cb_dtf_B) <- c("extinct", 
                        "pop_drop20", "pop_drop50", "pop_drop80",
                        "r_allele_peak15", "r_allele_peak45", "r_allele_peak75",
                        "final_r_allele", "final_pop_size", "final_inf_prev",
                        "max_r_allele", "time_max_r_allele",
                        "max_inf_prev", "time_last_zero_inf",
                        "min_pop", "time_min_pop",
                        "first_20", "first_50", "first_80",
                        "last_20", "last_50", "last_80",
                        "tot_20", "tot_50", "tot_80",
                        "at_K95", "firstK95",
                        "r_ts_d0", "parm_number")

dc_dat_B <- readRDS("dat/dc_fig_0911.Rdata")
dc_dtf_B <- as.data.frame(matrix(unlist(dc_dat_B), ncol = 29, byrow = T))
colnames(dc_dtf_B) <- c("extinct", 
                        "pop_drop20", "pop_drop50", "pop_drop80",
                        "r_allele_peak15", "r_allele_peak45", "r_allele_peak75",
                        "final_r_allele", "final_pop_size", "final_inf_prev",
                        "max_r_allele", "time_max_r_allele",
                        "max_inf_prev", "time_last_zero_inf",
                        "min_pop", "time_min_pop",
                        "first_20", "first_50", "first_80",
                        "last_20", "last_50", "last_80",
                        "tot_20", "tot_50", "tot_80",
                        "at_K95", "firstK95",
                        "r_ts_d0", "parm_number")

pm_dat_B <- readRDS("dat/pm_fig_0910.Rdata")
pm_dtf_B <- as.data.frame(matrix(unlist(pm_dat_B), ncol = 29, byrow = T))
colnames(pm_dtf_B) <- c("extinct", 
                        "pop_drop20", "pop_drop50", "pop_drop80",
                        "r_allele_peak15", "r_allele_peak45", "r_allele_peak75",
                        "final_r_allele", "final_pop_size", "final_inf_prev",
                        "max_r_allele", "time_max_r_allele",
                        "max_inf_prev", "time_last_zero_inf",
                        "min_pop", "time_min_pop",
                        "first_20", "first_50", "first_80",
                        "last_20", "last_50", "last_80",
                        "tot_20", "tot_50", "tot_80",
                        "at_K95", "firstK95",
                        "r_ts_d0", "parm_number")

#-----CALC PROBS & MERGE LIKE DATA: EX-----
str_ex_A <- str_dtf_A %>% filter(extinct == 1) %>% group_by(parm_number) %>% summarise(`P(extinct) - A` = n()/1000)
str_ex_B <- str_dtf_B %>% filter(extinct == 1) %>% group_by(parm_number) %>% summarise(`P(extinct) - B` = n()/1000)
str_ex <- merge(str_ex_A, str_ex_B, by = "parm_number", all = T)
str_ex$diff <- ifelse(is.na(str_ex$`P(extinct) - A`), 0, str_ex$`P(extinct) - A`) - ifelse(is.na(str_ex$`P(extinct) - B`), 0, str_ex$`P(extinct) - B`)

pm_ex_A <- pm_dtf_A %>% filter(extinct == 1) %>% group_by(parm_number) %>% summarise(`P(extinct) - A` = n()/1000)
pm_ex_B <- pm_dtf_B %>% filter(extinct == 1) %>% group_by(parm_number) %>% summarise(`P(extinct) - B` = n()/1000)
pm_ex <- merge(pm_ex_A, pm_ex_B, by = "parm_number", all = T)
pm_ex$diff <- ifelse(is.na(pm_ex$`P(extinct) - A`), 0, pm_ex$`P(extinct) - A`) - ifelse(is.na(pm_ex$`P(extinct) - B`), 0, pm_ex$`P(extinct) - B`)

dc_ex_A <- dc_dtf_A %>% filter(extinct == 1) %>% group_by(parm_number) %>% summarise(`P(extinct) - A` = n()/1000)
dc_ex_B <- dc_dtf_B %>% filter(extinct == 1) %>% group_by(parm_number) %>% summarise(`P(extinct) - B` = n()/1000)
dc_ex <- merge(dc_ex_A, dc_ex_B, by = "parm_number", all = T)
dc_ex$diff <- ifelse(is.na(dc_ex$`P(extinct) - A`), 0, dc_ex$`P(extinct) - A`) - ifelse(is.na(dc_ex$`P(extinct) - B`), 0, dc_ex$`P(extinct) - B`)

cb_ex_A <- cb_dtf_A %>% filter(extinct == 1) %>% group_by(parm_number) %>% summarise(`P(extinct) - A` = n()/1000)
cb_ex_B <- cb_dtf_B %>% filter(extinct == 1) %>% group_by(parm_number) %>% summarise(`P(extinct) - B` = n()/1000)
cb_ex <- merge(cb_ex_A, cb_ex_B, by = "parm_number", all = T)
cb_ex$diff <- ifelse(is.na(cb_ex$`P(extinct) - A`), 0, cb_ex$`P(extinct) - A`) - ifelse(is.na(cb_ex$`P(extinct) - B`), 0, cb_ex$`P(extinct) - B`)

fig_ex <- rbind(str_ex, pm_ex, dc_ex, cb_ex)

#-----CALC PROBS & MERGE LIKE DATA: ER-----
str_er_A <- str_dtf_A %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 &  at_K95 == 1 & r_allele_peak45 == 1) %>% group_by(parm_number) %>% summarise(`P(ER) - A` = n()/1000)
str_er_B <- str_dtf_B %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 &  at_K95 == 1 & r_allele_peak45 == 1) %>% group_by(parm_number) %>% summarise(`P(ER) - B` = n()/1000)
str_er <- merge(str_er_A, str_er_B, by = "parm_number", all = T)
str_er$diff <- ifelse(is.na(str_er$`P(ER) - A`), 0, str_er$`P(ER) - A`) - ifelse(is.na(str_er$`P(ER) - B`), 0, str_er$`P(ER) - B`)

pm_er_A <- pm_dtf_A %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 &  at_K95 == 1 & r_allele_peak45 == 1) %>% group_by(parm_number) %>% summarise(`P(ER) - A` = n()/1000)
pm_er_B <- pm_dtf_B %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 &  at_K95 == 1 & r_allele_peak45 == 1) %>% group_by(parm_number) %>% summarise(`P(ER) - B` = n()/1000)
pm_er <- merge(pm_er_A, pm_er_B, by = "parm_number", all = T)
pm_er$diff <- ifelse(is.na(pm_er$`P(ER) - A`), 0, pm_er$`P(ER) - A`) - ifelse(is.na(pm_er$`P(ER) - B`), 0, pm_er$`P(ER) - B`)

dc_er_A <- dc_dtf_A %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 &  at_K95 == 1 & r_allele_peak45 == 1) %>% group_by(parm_number) %>% summarise(`P(ER) - A` = n()/1000)
dc_er_B <- dc_dtf_B %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 &  at_K95 == 1 & r_allele_peak45 == 1) %>% group_by(parm_number) %>% summarise(`P(ER) - B` = n()/1000)
dc_er <- merge(dc_er_A, dc_er_B, by = "parm_number", all = T)
dc_er$diff <- ifelse(is.na(dc_er$`P(ER) - A`), 0, dc_er$`P(ER) - A`) - ifelse(is.na(dc_er$`P(ER) - B`), 0, dc_er$`P(ER) - B`)

cb_er_A <- cb_dtf_A %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 &  at_K95 == 1 & r_allele_peak45 == 1) %>% group_by(parm_number) %>% summarise(`P(ER) - A` = n()/1000)
cb_er_B <- cb_dtf_B %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 &  at_K95 == 1 & r_allele_peak45 == 1) %>% group_by(parm_number) %>% summarise(`P(ER) - B` = n()/1000)
cb_er <- merge(cb_er_A, cb_er_B, by = "parm_number", all = T)
cb_er$diff <- ifelse(is.na(cb_er$`P(ER) - A`), 0, cb_er$`P(ER) - A`) - ifelse(is.na(cb_er$`P(ER) - B`), 0, cb_er$`P(ER) - B`)

fig_er <- rbind(str_er, pm_er, dc_er, cb_er)

#-----CALC PROBS & MERGE LIKE DATA: IL-----
str_il_A <- str_dtf_A %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 0 & final_inf_prev < 0.25) %>% group_by(parm_number) %>% summarise(`P(IL) - A` = n()/1000)
str_il_B <- str_dtf_B %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 0 & final_inf_prev < 0.25) %>% group_by(parm_number) %>% summarise(`P(IL) - B` = n()/1000)
str_il <- merge(str_il_A, str_il_B, by = "parm_number", all = T)
str_il$diff <- ifelse(is.na(str_il$`P(IL) - A`), 0, str_il$`P(IL) - A`) - ifelse(is.na(str_il$`P(IL) - B`), 0, str_il$`P(IL) - B`)

pm_il_A <- pm_dtf_A %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 0 & final_inf_prev < 0.25) %>% group_by(parm_number) %>% summarise(`P(IL) - A` = n()/1000)
pm_il_B <- pm_dtf_B %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 0 & final_inf_prev < 0.25) %>% group_by(parm_number) %>% summarise(`P(IL) - B` = n()/1000)
pm_il <- merge(pm_il_A, pm_il_B, by = "parm_number", all = T)
pm_il$diff <- ifelse(is.na(pm_il$`P(IL) - A`), 0, pm_il$`P(IL) - A`) - ifelse(is.na(pm_il$`P(IL) - B`), 0, pm_il$`P(IL) - B`)

dc_il_A <- dc_dtf_A %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 0 & final_inf_prev < 0.25) %>% group_by(parm_number) %>% summarise(`P(IL) - A` = n()/1000)
dc_il_B <- dc_dtf_B %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 0 & final_inf_prev < 0.25) %>% group_by(parm_number) %>% summarise(`P(IL) - B` = n()/1000)
dc_il <- merge(dc_il_A, dc_il_B, by = "parm_number", all = T)
dc_il$diff <- ifelse(is.na(dc_il$`P(IL) - A`), 0, dc_il$`P(IL) - A`) - ifelse(is.na(dc_il$`P(IL) - B`), 0, dc_il$`P(IL) - B`)

cb_il_A <- cb_dtf_A %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 0 & final_inf_prev < 0.25) %>% group_by(parm_number) %>% summarise(`P(IL) - A` = n()/1000)
cb_il_B <- cb_dtf_B %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 0 & final_inf_prev < 0.25) %>% group_by(parm_number) %>% summarise(`P(IL) - B` = n()/1000)
cb_il <- merge(cb_il_A, cb_il_B, by = "parm_number", all = T)
cb_il$diff <- ifelse(is.na(cb_il$`P(IL) - A`), 0, cb_il$`P(IL) - A`) - ifelse(is.na(cb_il$`P(IL) - B`), 0, cb_il$`P(IL) - B`)

fig_il <- rbind(str_il, pm_il, dc_il, cb_il)

#-----PLOT HISTS-----
# add column and pivot longer
fig_ex$metric <- rep("P(extinct)")
fig_er$metric <- rep("P(ER)")
fig_il$metric <- rep("P(IL)")
fig_tot <- rbind(fig_ex[, c(1, 4, 5)], fig_er[, c(1, 4, 5)], fig_il[, c(1, 4, 5)])
fig_long <- pivot_longer(fig_tot, cols = 2) 

# box plots?
comp_sim_box <- ggplot(fig_long, aes(x = value, fill = metric)) + 
  geom_boxplot(col = "black") + 
  scale_fill_manual(values = c("#ac1457", "#DB6341", "#f1c4a2")) + 
  facet_wrap(~metric, nrow = 3) + 
  labs(x = "difference across two test simulations") + 
  theme_bw() + theme(legend.position = "none")

# histograms?
comp_sim_hist <- ggplot(fig_long, aes(x = value, fill = metric)) + 
  geom_histogram(col = "black", bins = 30) + 
  scale_fill_manual(values = c("#ac1457", "#DB6341", "#f1c4a2")) + 
  facet_wrap(~metric, nrow = 3) + 
  labs(x = "difference across two test simulations", y = NULL) + 
  theme_bw() + theme(legend.position = "none", 
                     axis.text.x = element_text(angle = 45, hjust = 1))

# save figure
png("figs/figure_plot/comp_sim_hist.png",height=170,width=85,res=400,units='mm')
print(comp_sim_hist)
dev.off()

quantile(fig_ex$diff, c(0.025, 0.975))
quantile(fig_er$diff, c(0.025, 0.975))
quantile(fig_il$diff, c(0.025, 0.975))
# use largest for the intervals