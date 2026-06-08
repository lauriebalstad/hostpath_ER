library(dplyr)
library(tidyverse)
# library(randomForest)
library(ggplot2)
# library(rpart)
# library(rpart.plot)
# library(optRF)
# library(ranger)
library(cowplot)
library(randomForestSRC)
# library(sensobol)

# loading the data -
sum_a <- readRDS("dat/gsa_result_0605-1547.Rdata")
sum_b <- readRDS("dat/gsa_result_0606-1626.Rdata")
sum_c <- readRDS("dat/gsa_result_0606-0502.Rdata")
sum_d <- readRDS("dat/gsa_result_0607-1214.Rdata")
sum_e <- readRDS("dat/gsa_result_0606-1934.Rdata")
random_parms <- readRDS("dat/mat_var_0604.Rdata") # figure out what these values were converging around....? new values have super little er
bd_parms <- readRDS("dat/mat_bd_0604.Rdata") # figure out what these values were converging around....? new values have super little er

sum_dat <- matrix(c(unlist(sum_a), 
                    unlist(sum_b), 
                    unlist(sum_c), 
                    unlist(sum_d), 
                    unlist(sum_e)), ncol = 27, byrow = T)
colnames(sum_dat) <- c("extinct", 
                       "pop_drop20", "pop_drop50", "pop_drop80",
                       "r_allele_peak15", "r_allele_peak45", "r_allele_peak75",
                       "final_r_allele", 
                       # "final_pop_size", 
                       "final_inf_prev",
                       "max_r_allele", # "time_max_r_allele",
                       "max_inf_prev", "time_last_zero_inf",
                       "min_pop", "time_min_pop",
                       "first_20", "first_50", "first_80",
                       "last_20", "last_50", "last_80",
                       "tot_20", "tot_50", "tot_80",
                       "at_K95", "firstK95",
                       "r_ts_d0", 
                       "parm_number")

# change -1s back to 0 bc that's their true rate for random_parms
random_parms <- as.data.frame(random_parms)
for (i in 1:dim(random_parms)[1]) {
  if (all(random_parms[i,3:5] == -2)) random_parms[i,3:5] <- 0 
}

# trying as a classifier
# set up data: merge parameters (V35) and simulated data (parm_number)
RF_class <- merge(random_parms, sum_dat, by.x = "V34", by.y = "parm_number")
RF_class <- merge(RF_class, bd_parms, by.x = "V34", by.y = "number")
# rename things (note need to keep these nice for RF)
colnames(RF_class) <- c("param_num", "compartments", "order", 
                        "B_F_RR", "B_F_WR", "B_F_WW", "SD_B_F", 
                        "B_D_RR", "B_D_WR", "B_D_WW", "SD_B_D", 
                        "M_S_RR", "M_S_WR", "M_S_WW", "M_I_RR", "M_I_WR", "M_I_WW", "SD_M_S", "SD_M_I", 
                        "G_RR", "G_WR", "G_WW", "SD_G", 
                        "L_RR", "L_WR", "L_WW", "mutation",
                        "K", "SD_K",
                        "init_time", "disease_gen", 
                        "init_allele", "init_inif", "final_time", "adaptive_cat", colnames(RF_class[36:65]))
# note the the variables that are categorical
RF_class$"compartments" <- as.factor(RF_class$"compartments")
RF_class$"order" <- as.factor(RF_class$"order")
RF_class$"adaptive_cat" <- as.factor(RF_class$"adaptive_cat")

# figure out classifications
RF_class$ER <- ifelse(RF_class$extinct == 0 & 
                        RF_class$pop_drop50 == 1 & 
                        abs(RF_class$last_50-RF_class$first_50-RF_class$tot_50) < 0.8*RF_class$tot_50 & 
                        RF_class$tot_50 > 3 & 
                        RF_class$at_K95 == 1 & 
                        RF_class$r_allele_peak45 == 1, 
                      1, 0) 

RF_class$T_ER <- ifelse(RF_class$extinct == 0 & 
                          RF_class$pop_drop50 == 1 & 
                          abs(RF_class$last_50-RF_class$first_50-RF_class$tot_50) < 0.8*RF_class$tot_50 & 
                          RF_class$tot_50 > 3 & 
                          RF_class$at_K95 == 1 & 
                          RF_class$r_allele_peak45 == 1 & 
                          RF_class$final_r_allele < 0.1 & RF_class$final_inf_prev < 0.1, 
                        1, 0) # <1% of outcomes

RF_class$IL <- ifelse(RF_class$extinct == 0 & 
                        RF_class$pop_drop50 == 1 & 
                        abs(RF_class$last_50-RF_class$first_50-RF_class$tot_50) < 0.8*RF_class$tot_50 & 
                        RF_class$tot_50 > 3 & 
                        RF_class$at_K95 == 1 & 
                        RF_class$r_allele_peak45 == 0 & 
                        RF_class$final_inf_prev == 0, 
                      1, 0) 

# have extinction column already
RF_class$Ext <- RF_class$extinct 

RF_class$NR <- ifelse(RF_class$extinct == 0 &
                      RF_class$pop_drop50 == 0 & 
                      RF_class$at_K == 1, 1, 0) # is very small!

RF_class$clss <- ifelse(RF_class$Ext == 1, "Ext", 
                        ifelse(RF_class$ER == 1, "ER", 
                               ifelse(RF_class$IL == 1, "IL",
                                      ifelse(RF_class$NR == 1, "NR", "UNK")))) # about 9% unknown

RF_class$clss <- as.factor(RF_class$clss) 
# plot(RF_class$clss) # mostly ER, extinction, some NR/UNK

# try distributed random forest, which deals with the fact that the outcome is a distribution not a category
# make outcome matrix
# get probabilities
p_ex <- RF_class %>% group_by(param_num) %>% summarise(p_ex = sum(Ext)/2500)
p_er <- RF_class %>% group_by(param_num) %>% summarise(p_er = sum(ER)/2500)
p_il <- RF_class %>% group_by(param_num) %>% summarise(p_il = sum(IL)/2500)
p_nr <- RF_class %>% group_by(param_num) %>% summarise(p_nr = sum(NR)/2500)
p_uk <- RF_class %>% group_by(param_num) %>% summarise(p_uk = sum(clss=="UNK")/2500) 
# merge everyone
tmp <- merge(p_ex, p_er, all = T)
tmp <- merge(tmp, p_il, all = T)
tmp <- merge(tmp, p_nr, all = T)
out_mat <- merge(tmp, p_uk, all = T)
# check that out_mat is adding to 1
tmp_check <- out_mat %>% mutate(tot_prob = p_ex + p_er + p_il + p_nr + p_uk); range(tmp_check$tot_prob)
# clean up the input data -- only need one row for each now
pred_vars <- c(2, 3, 6, 10, 14, 17, 22, 26:35, 62:65)
in_mat <- RF_class[seq(from = 1, to = dim(RF_class)[1], by = 2500), pred_vars]

# now run the forest
tmp_df <- cbind(in_mat, out_mat[,2:6]) 
clss_mv <- rfsrc(Multivar(p_ex, p_er, p_il, p_nr, p_uk) ~ . + compartments*adaptive_pathway, data = tmp_df, # 
                 importance = "permute", ntree = 500)
plot(clss_mv) # first look for one outcome only!
imp <- vimp(clss_mv, importance ="permute") # this repeats variable calc from above fyi
# use ``imp$regrOutput$p_ex$importance`` to get importance for a particular outcome
# plot.variable(clss_mv, m.target = "p_er") 
saveRDS(clss_mv, "sens/rf_output.Rdata")

# try predicting a few rows of fig1B, to check out how off base things are
# merge data w/original cases
str_dat <- readRDS("dat/str_fig_0605.Rdata")
# convert to df
str_dtf <- as.data.frame(matrix(unlist(str_dat), ncol = 27, byrow = T))
colnames(str_dtf) <- c("extinct", 
                       "pop_drop20", "pop_drop50", "pop_drop80",
                       "r_allele_peak15", "r_allele_peak45", "r_allele_peak75",
                       "final_r_allele", 
                       # "final_pop_size", 
                       "final_inf_prev",
                       "max_r_allele", # "time_max_r_allele",
                       "max_inf_prev", "time_last_zero_inf",
                       "min_pop", "time_min_pop",
                       "first_20", "first_50", "first_80",
                       "last_20", "last_50", "last_80",
                       "tot_20", "tot_50", "tot_80",
                       "at_K95", "firstK95",
                       "r_ts_d0", 
                       "parm_number")

# merge data w/original cases
cases <- expand.grid("transmission" = c("Density only", "Density w/reservior", "Lasting disease shock"), 
                     "compartments" = c("Mortality", "Recovery", "Immunity"), # note SIR doesn't seem to have a big effect
                     "robustness" = c("MB", "TB", "RA", "N")) # 1 = mortality, 2 = transmission, 3 = recovery, 4 = demographic rescue only
# remove 1/3 compartment/robustness combo
cases <- cases[c(1:18, 22:36), ]
N <- dim(cases)[1]
cases$number <- 1:N
# get summaries of each probability to plot
prob_dat <- str_dtf %>% group_by(parm_number) %>% 
  summarize(`P(Ext)` = sum(extinct)/n()) 
er_dat <- str_dtf %>% 
  group_by(parm_number) %>% 
  filter(extinct == 0 & abs(last_50-first_50-tot_50+1) < 0.8*tot_50 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 1) %>% 
  summarize(`P(All ER)` = n()/2500)
il_dat <- str_dtf %>% filter(extinct == 0 & abs(last_50-first_50-tot_50+1) < 0.8*tot_50 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 0 & final_inf_prev == 0) %>% 
  group_by(parm_number) %>% 
  summarize(`P(IL)` = n()/2500)
# add in a no risk category
nr_dat <- str_dtf %>% filter(extinct == 0 & pop_drop50 == 0 & at_K95 == 1) %>% 
  group_by(parm_number) %>% 
  summarize(`P(NR)` = n()/2500)
# merge
tmp <- merge(prob_dat, er_dat, by = c("parm_number"), all = T)
tmp <- merge(tmp, il_dat, by = c("parm_number"), all = T)
tmp <- merge(tmp, nr_dat, by = c("parm_number"), all = T)
plot_dat <- merge(tmp, cases, by.x = c("parm_number"), by.y = "number")
plot_dat[is.na(plot_dat)] <- 0 # replace things

# now use predict to predict what those probabilities would be from the RF
base_vect <- c(2, 4, # compartments, event order
               -1, -1, # freq & density transmission rate WW
               0.05, 1, # mortality WW
               0.05,  # (7) recovery WW 
               1, # reproduction WW 
               0.005, # mutation
               100, 4, # (10) pop size w/sds
               120, 
               4, # disease gen
               0.05, 0.1, # allele & inf inits
               160, # final time
               0, 1, # pathway, benefit
               0.5, # dominance
               0.05, # (20) general sd
               0.7) # cost
pred_results_str <- NULL
rep_cases <- cases %>% filter(transmission != "Lasting disease shock")
N <- dim(rep_cases)[1]
# loop through w predict
for (i in 1:N) {
  
  case_vect <- base_vect
  
  # set disease transmission
  if (rep_cases$transmission[i]=="Density only") {
    case_vect[4] <- 1 # transmission increased (dens)
    case_vect[3] <- 0 # by def for this figure/RF
  }
  if (rep_cases$transmission[i]=="Density w/reservior") {
    case_vect[3] <- 0.05 # transmission increased (freq)
    case_vect[4] <- 1
  }
  if (rep_cases$robustness[i] == 1) {
    case_vect[17] <- 1
    case_vect[18] <- 0 # benefit is defined as the wild type*benefit for robust type (100% reduction here)
  } # 1 = mortality
  if (rep_cases$robustness[i] == 2) {
    case_vect[17] <- 2
    case_vect[18] <- 0 # benefit is defined as the wild type*benefit for robust type (100% reduction here)
  } # 1 = transmission
  if (rep_cases$robustness[i] == 3) {
    case_vect[17] <- 3
    case_vect[18] <- 0.5 # benefit is defined as the wild type*benefit for robust type (100% reduction here)
  } # 1 = recovery, nb: recovery higher
  
  # modify base parameters -- compartment
  case_vect[1] <- rep_cases$compartments[i]
  
  case_vect <- as.data.frame(t(case_vect)); colnames(case_vect) <- clss_mv[["xvar.names"]]
  
  tmp <- predict(clss_mv, case_vect)
  
  out_dat <- c(tmp$regrOutput$p_er$predicted,
               tmp$regrOutput$p_ter$predicted,
               tmp$regrOutput$p_ex$predicted, 
               tmp$regrOutput$p_il$predicted, 
               tmp$regrOutput$p_nr$predicted, 
               tmp$regrOutput$p_uk$predicted
  )
  
  pred_results_str <- rbind(pred_results_str, t(out_dat))
  
}
# organize so it's matching w plot_long
pred_results_str <- as.data.frame(pred_results_str)
pred_results_str <- cbind(rep_cases$number, pred_results_str)
colnames(pred_results_str) <- c("parm_number", "P(All ER)", "P(Ext)", "P(IL)", "P(NR)", "Inconclusive")
# double check that things are adding to 1
pred_check <- apply(pred_results_str[, 2:6], 1, sum); range(pred_check)
pred_dat <- merge(pred_results_str, cases, by.x = c("parm_number"), by.y = "number")

# pivot longer for both
pred_dat$shp <- rep("RF predicted")
plot_dat$shp <- rep("full sim.")
plot_long <- pivot_longer(plot_dat, cols = 2:5, names_to = "outcome", values_to = "value")
pred_long <- pivot_longer(pred_dat, cols = 2:6, names_to = "outcome", values_to = "value") # ignore unks
dat_long <- rbind(plot_long, pred_long)

# points plot -- note RF doesn't know what N is, so drop that from this comp
comp_sim_pred <- ggplot(dat_long %>% filter(as.character(robustness) != "N" & 
                                            transmission != "Lasting disease shock"), 
                        aes(robustness, value, col = shp, pch = transmission)) + 
  geom_point(size = 2.25, position = position_dodge(width = 0.7)) + 
  scale_color_manual(values = c("black", "gray80")) +
  scale_shape_manual(values = c(15, 16)) + 
  facet_grid(rows = vars(outcome), cols = vars(compartments)) + 
  labs(x = "Adaptive pathway", y = "Probability", col = "Calc. source", pch = "Transmission") + 
  theme_bw(base_size = 12) + 
  theme(legend.position = "bottom", 
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.key.spacing.y = unit(0, "pt")) + 
  labs(x = "Mutant allele benefit\n(percent change in rate)", y = NULL)

pdf("figs/figure_plot/figS4_HP.pdf",height=175/25.4,width=170/25.4)
print(comp_sim_pred)
dev.off()

# using information from 02_run_randomForest
RF_class$clss <- recode(RF_class$clss, "IL" = "Inf. Loss", "UNK" = "Inconclusive")
figS2 <- ggplot(RF_class, aes(x = clss, fill = clss)) + geom_bar() + 
  scale_fill_manual(values = c("#DB6341", "#ac1457", "#f1c4a2", "gray40", "gray80")) +
  theme_bw(base_size = 12) + labs(x = "Classification") + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

pdf("figs/figure_plot/figS2_HP.pdf",height=90/25.4,width=85/25.4)
print(figS2)
dev.off()

