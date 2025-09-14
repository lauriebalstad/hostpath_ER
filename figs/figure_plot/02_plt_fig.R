library(tidyverse)
library(cowplot)
library(ggplot2)
library(viridis)
library(ggh4x)

#-----READ IN DATA-----
# str_dat <- readRDS("C:/Users/lauri/AppData/Roaming/MobaXterm/home/evorescue_hostpath_str/dat/str_fig_d2_0203.Rdata")
# str_dat <- readRDS("dat/str_fig_d2_0203.Rdata")
str_dat <- readRDS("dat/str_fig_0909.Rdata")
# convert to df
str_dtf <- as.data.frame(matrix(unlist(str_dat), ncol = 29, byrow = T))
colnames(str_dtf) <- c("extinct", 
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

# cb_dat <- readRDS("C:/Users/lauri/AppData/Roaming/MobaXterm/home/evorescue_hostpath_str/dat/cb_fig_d2_0203.Rdata")
# cb_dat <- readRDS("dat/cb_fig_d2_0203.Rdata")
cb_dat <- readRDS("dat/cb_fig_0909.Rdata")
cb_dtf <- as.data.frame(matrix(unlist(cb_dat), ncol = 29, byrow = T))
colnames(cb_dtf) <- c("extinct", 
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

# dc_dat <- readRDS("C:/Users/lauri/AppData/Roaming/MobaXterm/home/evorescue_hostpath_str/dat/dc_fig_d2_0203.Rdata")
# dc_dat <- readRDS("dat/dc_fig_d2_0203.Rdata")
dc_dat <- readRDS("dat/dc_fig_0910.Rdata")
dc_dtf <- as.data.frame(matrix(unlist(dc_dat), ncol = 29, byrow = T))
colnames(dc_dtf) <- c("extinct", 
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

# 0205 looks nice as well pm_fig_d2_0203
# pm_dat <- readRDS("C:/Users/lauri/AppData/Roaming/MobaXterm/home/evorescue_hostpath_str/dat/pm_fig_d2_0203.Rdata")
# pm_dat <- readRDS("dat/pm_fig_d2_0203.Rdata")
pm_dat <- readRDS("dat/pm_fig_0909.Rdata")
pm_dtf <- as.data.frame(matrix(unlist(pm_dat), ncol = 29, byrow = T))
colnames(pm_dtf) <- c("extinct", 
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

#-----FIG1: TIME SERIES & STR-----
# time series data loading
str_dat <- readRDS("dat/all_ts_dat_0204.rds") # conintue to play w... there's a lot of cases w/o pop drop getting hit
# rotate data
# remove nans as 0s
# str_dat$`Infection\nprevelence` <- ifelse(is.nan(str_dat$`Infection\nprevelence`), 0, str_dat$`Infection\nprevelence`)
# remove exs without popdrop
kps <- unique((str_dat %>% filter(`A - Scaled\npopulation size` < 0.5))$trial)
all_ts_long <- pivot_longer(str_dat %>% filter(trial %in% kps) %>% select(c(1, 7, 9, 11)), cols = c(3:5), names_to = "outcome", values_to = "value")
all_ts_long$outcome_f <- factor(all_ts_long$outcome, 
                                levels = c("A - Scaled\npopulation size", "B - Adaptive allele\nfrequency", "C - Infection\nprevelence"), 
                                labels = c("Scaled pop. size", "Adaptive allele freq.", "Infection prev."))
colors <- c("Extinction" = alpha("black", 0.8), "Persist. ER" = alpha("#ac1457",0.8),
            "Temp. ER" = "#f1c4a2", "Lose Inf." = alpha("#DB6341", 0.8))
fig1A <- ggplot(all_ts_long, aes(ts, value)) + # 7, 8, 9, 10 gives exinction
  geom_line(aes(group = trial), col = alpha("black", 0.1), lwd = 0.75) +
  geom_line(data = all_ts_long %>% filter(trial == 6), aes(color = "Extinction"), lwd = 0.85) + # extinction
  geom_line(data = all_ts_long %>% filter(trial == 9), aes(color = "Persist. ER"), lwd = 0.85) + # ER
  geom_line(data = all_ts_long %>% filter(trial == 20), aes(color = "Lose Inf."), lwd = 0.85) + # DR
  geom_line(data = all_ts_long %>% filter(trial == 19), aes(color = "Temp. ER"), lwd = 0.85) + # back selection
  facet_wrap(~outcome_f, scale = "free_y", ncol = 1) + theme_bw() + 
  scale_color_manual(values = colors) + 
  theme(text = element_text(size = 12), legend.position = "bottom") + 
  guides(color = guide_legend(nrow = 2)) +
  labs(x = "Standardized generation", y = NULL, col = NULL)

# merge data w/original cases
cases <- expand.grid("transmission type" = 1:2, 
                     "compartments" = 1:3, # note SIR doesn't seem to have a big effect
                     "robustness" = 1:4) # 1 = mortality, 2 = transmission, 3 = recovery, 4 = demographic rescue only
# remove 1/3 compartment/robustness combo
cases <- cases[c(1:12, 15:24), ]
N <- dim(cases)[1]
cases$number <- 1:N
tmp <- merge(str_dtf, cases, by.x = "parm_number", by.y = "number", all = T)
# get summaries to plot
prob_dat <- tmp %>% group_by(parm_number) %>% 
  summarize(`P(Ext)` = sum(extinct)/n()) 
# er_dat <- tmp %>% filter(extinct == 0 & pop_drop50 == 1 & at_K95 == 1 & r_allele_peak45 == 1) %>% 
#   group_by(parm_number) %>% 
#   summarize(`P(ER)` = n()/1000)
er_dat <- tmp %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 1) %>% 
  group_by(parm_number) %>% 
  summarize(`P(ER)` = n()/1000)
il_dat <- tmp %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 0 & final_inf_prev < 0.25) %>% 
  group_by(parm_number) %>% 
  summarize(`P(IL)` = n()/1000)
# merge
tmp <- merge(prob_dat, er_dat, by = c("parm_number"), all = T)
tmp <- merge(tmp, il_dat, by = c("parm_number"), all = T)
plot_dat <- merge(tmp, cases, by.x = c("parm_number"), by.y = "number")
# rename
plot_dat$compartments <- recode(plot_dat$compartments, "1" = "SIX", "2" = "SIS", "3" = "SIR")
plot_dat$`transmission type` <- recode(plot_dat$`transmission type`, "1" = "density", "2" = "environmental", "3" = "density + environmental")
plot_dat$robustness <- recode(plot_dat$robustness, "1" = "\u03bc", "2" = "\u03b2", "3" = "\u03d2", "4" = "N")
# convert to long
plot_long <- pivot_longer(plot_dat, cols = 2:4, names_to = "outcome", values_to = "value")
plot_long$value <- ifelse(is.na(plot_long$value), 0, plot_long$value)
plot_long$plus <- plot_long$value + 0.036
plot_long$less <- plot_long$value - 0.036
for (i in 1: length(plot_long$value)) {
  if (plot_long$plus[i] > 1) plot_long$plus[i] <- 1
  if (plot_long$less[i] < 0) plot_long$less[i] <- 0
}
plot_long$outcome_f <- factor(plot_long$outcome, levels = c("P(Ext)", "P(ER)", "P(IL)"))
# points plot
fig1B <- ggplot(plot_long, aes(robustness, value, col = `transmission type`)) + 
  geom_linerange(aes(ymin = less, ymax = plus), lwd = 1) + 
  geom_hline(data = plot_long %>% filter(robustness == "N"), aes(yintercept = value), col = "gray70", lty = "dashed") + 
  geom_point(size = 2) + scale_color_manual(values = c("#ac1457", "#f1c4a2")) +
  facet_grid(rows = vars(outcome_f), cols = vars(compartments)) + 
  labs(x = "Robustness type", y = "Probability", col = "Transmission type") + 
  theme_bw()  + 
  theme(text = element_text(size = 12), legend.position = "bottom") 

# combine figure 1
fig1 <- plot_grid(fig1A, fig1B, rel_widths= c(0.65, 1), nrow = 1)

# save figure
png("figs/figure_plot/structural_overview_0909B.png",height=145,width=170,res=400,units='mm')
print(fig1)
dev.off()

#-----FIG2: COST/BENEFIT-----
# merge data
sis_cb <- expand.grid("benefit" = c(0, 0.25, 0.5, 1), # 0 = total benefit, 1.5 = some cost for B & M, but opposite for G
                      "cost" = c(1.7, 1.9, 2.1), # 1.7 = some cost, 2.1 = no cost
                      # "cost" = c(1.1, 1.9, 2.1), # 1 = some cost, 2 = no cost
                      "case" = c("\u03b2", "\u03bc", "\u03d2"), 
                      "compartments" = c(2, 3)) 
N <- dim(sis_cb)[1]
sis_cb$number <- 1:N
cost_ben <- merge(cb_dtf, sis_cb, by.x = "parm_number", by.y = "number", all = T)
# cost_ben <- cost_ben %>% filter(cost %in% c(1, 2)) # cost == 1 is okay to plot too, just adds unneeded points
cost_ben <- cost_ben %>% filter(compartments == 2)
# get summaries to plot -- P(ER) + T_K + max/final R/inf (6) for the ER cases only
ex_dat <- cost_ben %>% filter(extinct == 1) %>% 
  group_by(benefit, cost, case, compartments) %>% 
  summarize(`P(Ext)` = n()/1000)
er_dat <- cost_ben %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 1) %>% 
  group_by(benefit, cost, case, compartments) %>% 
  summarize(`P(ER)` = n()/1000)
dr_dat <- cost_ben %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 0) %>%    
  group_by(benefit, cost, case, compartments) %>% 
  summarize(`P(IL)` = n()/1000)
plt_dat <- merge(ex_dat, er_dat) # , dr_dat)
# rbind creates a ton of NAs but it's okay...
cb_summ_long <- pivot_longer(plt_dat, cols= 5:6, names_to = "outcome", values_to = "value")
cb_summ_long$plus <- cb_summ_long$value + 0.036
cb_summ_long$less <- cb_summ_long$value - 0.036
for (i in 1:length(cb_summ_long$value)) {
  if (cb_summ_long$plus[i] > 1) cb_summ_long$plus[i] <- 1
  if (cb_summ_long$less[i] < 0) cb_summ_long$less[i] <- 0
}
cost_ben$firstK95 <- cost_ben$firstK95-45 # remove 100 year burn in period (100/2 disease gens)
cost_ben <- rename(cost_ben, 
                 `Final R allele\nfrequency` = final_r_allele, 
                 `Final infection\nprevelence` = final_inf_prev, 
                 `Time to\nK (ER)` = firstK95)
cb_all_long <- pivot_longer(cost_ben %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 1) %>% 
                              select(c(9, 11, 28, 30:33)), 
                            cols = 1:3, names_to = "outcome", values_to = "value")
cb_summ_long$outcome_f <- factor(cb_summ_long$outcome, levels = c("P(Ext)", "P(ER)", "P(IL)", 
                                                            "Final R allele\nfrequency", "Final infection\nprevelence", "Time to\nK (ER)"))
cb_all_long$outcome_f <- factor(cb_all_long$outcome, levels = c("P(Ext)", "P(ER)", "P(IL)", 
                                                            "Final R allele\nfrequency", "Final infection\nprevelence", "Time to\nK (ER)"))
# remove the no-evo cases from box plots
cb_mu_beta <- cb_all_long %>% filter(case %in% c("\u03bc", "\u03b2")) %>% filter(benefit < 1)
cb_gamma <- cb_all_long %>% filter(case == "\u03d2") %>% filter(benefit > 0)
cb_no_evo <- rbind(cb_mu_beta, cb_gamma)
# only 1.7/2.1
cb_no_evo_HL <- cb_no_evo %>% filter(cost %in% c(1.7, 2.1))
cb_summ_long_HL <- cb_summ_long %>% filter(cost %in% c(1.7, 2.1))
fig2 <- ggplot(NULL, aes(x=as.factor(benefit), y=value, col = as.factor(cost))) + 
  geom_boxplot(data = cb_no_evo_HL) + 
  geom_linerange(data = cb_summ_long_HL, aes(ymin = less, ymax = plus), lwd = 1) + 
  geom_point(data = cb_summ_long_HL, size = 2) + 
  scale_color_manual(values = c("#ac1457", "#f1c4a2")) + 
  facet_grid(rows = vars(outcome_f), cols = vars(case), scales = "free_y") + 
  theme_bw() + theme(text = element_text(size = 12), legend.position = "bottom") + 
  labs(x = "Allele strength", col = "Allele fecundity", y = NULL)
# save figure
png("figs/figure_plot/cost_benefit_0909B.png",height=175,width=170,res=400,units='mm')
print(fig2)
dev.off()

#-----FIG3: POP SIZES-----
# merge data
sis_pop_mut <- expand.grid("population size" = c(50, 100, 500), # probably provides enough insight...
                           "robustness" = c(1, 3), # M, G only bc they are strongest ER (from fig 3)
                           "mutation rate" = c(0.0005, 0.005, 0.01),
                           "compartments" = c(2, 3)) # just for variation
N <- dim(sis_pop_mut)[1]
sis_pop_mut$number <- 1:N
ex_risk <- merge(pm_dtf, sis_pop_mut, by.x = "parm_number", by.y = "number", all = T)
ex_risk <- ex_risk %>% filter(compartments == 2, robustness == 3)
# rename
ex_risk$robustness <- recode(ex_risk$robustness, "1" = "\u03bc", "2" = "\u03b2", "3" = "\u03d2", "4" = "N")
ex_risk$compartments <- recode(ex_risk$compartments, "1" = "SIX", "2" = "SIS", "3" = "SIR")
# get summaries to plot -- P(ER) + T_K + max/final R/inf (6) for the ER cases only
ex_dat <- ex_risk %>% filter(extinct == 1) %>% 
  group_by(parm_number) %>% 
  summarize(`P(Ext)` = n()/1000)
er_dat <- ex_risk %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 1) %>% 
  # group_by(`population size`, `robustness`, `mutation rate`, `compartments`) %>% 
  group_by(parm_number) %>% 
  summarize(`P(ER))` = n()/1000)
plt_dat <- merge(ex_dat, er_dat, by = "parm_number", all = T)
tmp <- merge(plt_dat, sis_pop_mut, by.x = "parm_number", by.y = "number")
tmp$compartments <- recode(tmp$compartments, "1" = "SIX", "2" = "SIS", "3" = "SIR")
tmp$robustness <- recode(tmp$robustness, "1" = "\u03bc", "2" = "\u03b2", "3" = "\u03d2", "4" = "N")
# remove NAs --> those are 0s
tmp[is.na(tmp)] <- 0
# plt_dat <- rbind(ex_dat, er_dat)
pm_summ_long <- pivot_longer(tmp[, 2:7], cols= 1:2, names_to = "outcome", values_to = "value")
pm_summ_long$plus <- pm_summ_long$value + 0.036
pm_summ_long$less <- pm_summ_long$value - 0.036
for (i in 1:length(pm_summ_long$value)) {
  if (pm_summ_long$plus[i] > 1) pm_summ_long$plus[i] <- 1
  if (pm_summ_long$less[i] < 0) pm_summ_long$less[i] <- 0
}
ex_risk$firstK95 <- ex_risk$firstK95-45 # remove 50 year burn in period
ex_risk <- rename(ex_risk, 
                   `Time to\nK (ER)` = firstK95)
pm_all_long <- pivot_longer(ex_risk %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 1) %>% 
                              select(c(28, 30:33)), 
                            cols = 1, names_to = "outcome", values_to = "value")
custom_y <- list(
  scale_y_continuous(limits = c(0, 1)),
  scale_y_continuous(limits = c(0, 1)),
  scale_y_continuous(limits = c(0, 140))
)
pm_summ_long$outcome_f <- factor(pm_summ_long$outcome, levels = c("P(Ext)", "P(ER)", "P(IL)", 
                                                                  "Final R allele\nfrequency", "Final infection\nprevelence", "Time to\nK (ER)"))
pm_all_long$outcome_f <- factor(pm_all_long$outcome, levels = c("P(Ext)", "P(ER)", "P(IL)", 
                                                                "Final R allele\nfrequency", "Final infection\nprevelence", "Time to\nK (ER)"))
fig3 <- ggplot(NULL, aes(x=as.factor(`population size`), y=value, col=as.factor(`mutation rate`))) + 
  geom_boxplot(data = pm_all_long) + 
  # coord_cartesian(ylim = c(0, NA)) + 
  geom_linerange(data = pm_summ_long, aes(ymin = less, ymax = plus), lwd = 1) + 
  geom_point(data = pm_summ_long, size = 2) + 
  scale_color_manual(values = c("#ac1457","#DB6341", "#f1c4a2")) + 
  facet_grid(outcome ~ robustness, scales = "free_y") + 
  theme_bw() + theme(text = element_text(size = 12), legend.position = "bottom") + 
  labs(x = "Population size", col = "Mutation rate", y = NULL)
# save figure
png("figs/figure_plot/pop_mut_0909B.png",height=145,width=85,res=400,units='mm')
print(fig3)
dev.off()
# if 3x3 grid, width = 170 (this is pm_fig_A_0209 or pm_fig_B_0209, using mortality blocking (2))

#-----SUPP: DISEASE CYCLES-----
# load in data
sis_dc <- expand.grid("disease cycles" = 1:3, # probably provides enough insight...
                      "robustness" = c(1, 3), # M, G only bc they are strongest ER (from fig 3)
                      "event order" = c(2, 4, 6), 
                      "transmission" = c(1, 2),  # 1 = density, 2 = envrionemntal 
                      "compartments" = c(2, 3)) # 2 = SIS, 3 = SIR
N <- dim(sis_dc)[1]
sis_dc$parm_number <- 1:N
dc_dat <- merge(dc_dtf, sis_dc, by = "parm_number", all = T)
dc_dat <- dc_dat %>% filter(transmission == 1 & compartments == 2)
sis_dc <- sis_dc %>% filter(transmission == 1 & compartments == 2)
# dc_dat <- dc_dat %>% filter(robustness == 3 & compartments == 2)
# sis_dc <- sis_dc %>% filter(robustness == 3 & compartments == 2)
# rename
dc_dat$`event order` <- recode(dc_dat$`event order`, "2" = "\u03b2, \u03d2, \u03bc", "4" = "\u03bc, \u03b2, \u03d2", "6" = "\u03d2, \u03bc, \u03b2")
dc_dat$robustness <- recode(dc_dat$robustness, "1" = "\u03bc", "2" = "\u03b2", "3" = "\u03d2", "4" = "N")
dc_dat$compartments <- recode(dc_dat$compartments, "1" = "SIX", "2" = "SIS", "3" = "SIR")
dc_dat$transmission <- recode(dc_dat$transmission, "1" = "Density-dependent", "2" = "Environmental")
sis_dc$`event order` <- recode(sis_dc$`event order`, "2" = "\u03b2, \u03d2, \u03bc", "4" = "\u03bc, \u03b2, \u03d2", "6" = "\u03d2, \u03bc, \u03b2")
sis_dc$robustness <- recode(sis_dc$robustness, "1" = "\u03bc", "2" = "\u03b2", "3" = "\u03d2", "4" = "N")
sis_dc$compartments <- recode(sis_dc$compartments, "1" = "SIX", "2" = "SIS", "3" = "SIR")
sis_dc$transmission <- recode(sis_dc$transmission, "1" = "Density-dependent", "2" = "Environmental")
# get summaries to plot -- P(ER) + T_K + max/final R/inf (6) for the ER cases only --> note more generous def of pop drop
er_dat <- dc_dat %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 1) %>% 
  group_by(`parm_number`) %>% 
  summarize(`P(ER)` = n()/1000)
tmp <- merge(er_dat, sis_dc, by = "parm_number", all = T)
er_dat <- tmp; er_dat[is.na(er_dat)] <- 0
dr_dat <- dc_dat %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 0) %>%    
  group_by(parm_number) %>% 
  summarize(`P(IL)` = n()/1000)
tmp <- merge(dr_dat, sis_dc, by = "parm_number", all = T)
dr_dat <- tmp; dr_dat[is.na(dr_dat)] <- 0
ex_dat <- dc_dat %>% filter(extinct == 1) %>% 
  group_by(parm_number) %>% 
  summarize(`P(Ext)` = n()/1000)
tmp <- merge(ex_dat, sis_dc, by = "parm_number", all = T)
ex_dat <- tmp; ex_dat[is.na(ex_dat)] <- 0
# time to K
dc_dat$firstK95 <- dc_dat$firstK95-45
dc_dat <- rename(dc_dat, 
                   `Final R allele\nfrequency` = final_r_allele, 
                   `Final infection\nprevelence` = final_inf_prev, 
                   `Time to\nK (ER)` = firstK95)
dc_timeK_long <- pivot_longer(dc_dat %>% 
                                filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 1) %>% 
                                select(c(9, 11, 28, 30:34)), 
                              cols = 1:3, names_to = "outcome", values_to = "value")
tmp_dat <- merge(er_dat, ex_dat)
tmp_dat <- merge(tmp_dat, dr_dat)
# dc_er_long <- pivot_longer(er_dat, cols= c(2), names_to = "outcome", values_to = "value")
# dc_ex_long <- pivot_longer(ex_dat, cols= c(2), names_to = "outcome", values_to = "value")
# dc_dr_long <- pivot_longer(dr_dat, cols= c(2), names_to = "outcome", values_to = "value")
dc_all_long <- pivot_longer(tmp_dat[, 2:9], cols= 6:8, names_to = "outcome", values_to = "value")
dc_all_long$plus <- dc_all_long$value + 0.036
dc_all_long$less <- dc_all_long$value - 0.036
for (i in 1:length(dc_all_long$value)) {
  if (dc_all_long$plus[i] > 1) dc_all_long$plus[i] <- 1
  if (dc_all_long$less[i] < 0) dc_all_long$less[i] <- 0
}
dc_timeK_long$outcome_f <- factor(dc_timeK_long$outcome, levels = c("P(Ext)", "P(ER)", "P(IL)", 
                                                                  "Final R allele\nfrequency", "Final infection\nprevelence", "Time to\nK (ER)"))
dc_all_long$outcome_f <- factor(dc_all_long$outcome, levels = c("P(Ext)", "P(ER)", "P(IL)", 
                                                                "Final R allele\nfrequency", "Final infection\nprevelence", "Time to\nK (ER)"))
figS1 <- ggplot(NULL, aes(x=as.factor(`disease cycles`), y=value, 
                                  col=as.factor(robustness))) + 
  geom_boxplot(data = dc_timeK_long) + 
  geom_linerange(data = dc_all_long, aes(ymin = less, ymax = plus), lwd = 1) + 
  geom_point(data = dc_all_long, size = 2) + 
  scale_color_manual(values = c("#ac1457", "#f1c4a2")) + 
  facet_grid(row = vars(outcome_f), col = vars(`event order`), scales = "free_y") + 
  theme_bw() + theme(text = element_text(size = 12), legend.position = "bottom") + 
  labs(x = "Disease cycles between reproduction", col = "Robustness type", y = NULL)
# check labels/facet options! can either filter by one of the robustness or by one of the transmission 
# save figure
png("figs/figure_plot/dc_0909B_dens_sis.png",height=205,width=170,res=400,units='mm')
print(figS1)
dev.off()
