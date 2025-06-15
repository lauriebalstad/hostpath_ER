library(dplyr)
library(randomForest)
library(ggplot2)

# loading the data -- note weird directory because of HPC use
sum_0302a <- readRDS("dat/sim_dat_0302.Rdata")
sum_0302b <- readRDS("dat/sim_dat_2_0302.Rdata")
random_parms <- readRDS("dat/mat_var_0211.Rdata") # figure out what these values were converging around....? new values have super little er
# ^ read correct matrix!

# convert to df
sim_dat <- rbind(sum_0302a, sum_0302b) # just one???
# colnames(sim_dat) <- c("extinct", 
#                        "pop_drop20", "pop_drop50", "pop_drop80",
#                        "r_allele_peak15", "r_allele_peak45", "r_allele_peak75",
#                        "final_r_allele", "final_pop_size", "final_inf_prev",
#                        "max_r_allele", "time_max_r_allele",
#                        "max_inf_prev", "time_last_zero_inf",
#                        "min_pop", "time_min_pop",
#                        "at_K95", "firstK95",
#                        "r_ts_d0", "parm_number")
# summarize by case
rep_num <- 1000 # 500  
sens_dat_ext <- sim_dat %>% group_by(parm_number) %>% 
  summarise(`P(extinct)` = sum(extinct)/rep_num, 
            `P(pop_drop20)` = sum(pop_drop20)/rep_num,
            `P(pop_drop50)` = sum(pop_drop50)/rep_num,
            `P(pop_drop80)` = sum(pop_drop80)/rep_num, 
            avg_final_r = mean(final_r_allele),
            avg_max_r = mean(max_r_allele),
            avg_final_inf = mean(final_inf_prev), 
            abg_max_inf = mean(max_inf_prev), 
            time_K = mean(firstK95))
sens_dat_erE <- sim_dat %>% 
  filter(extinct == 0, pop_drop20 == 1, at_K95 == 1, r_allele_peak75 == 1) %>%
  group_by(parm_number) %>% 
  summarise(`P(ER_20_75)` = n()/rep_num)
tmp <- merge(sens_dat_ext, sens_dat_erE, by = "parm_number", all = T)
sens_dat_erM <- sim_dat %>% 
  filter(extinct == 0, pop_drop50 == 1, at_K95 == 1, r_allele_peak45 == 1) %>%
  group_by(parm_number) %>% 
  summarise(`P(ER_50_45)` = n()/rep_num, ER_K95 = mean(firstK95))
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

# remove NAs --> those are 0s
GSA_summary_dat[is.na(GSA_summary_dat)] <- 0

# get data together -- watch which variable merge is occuring w -- 35 for pre-0122, 36 for post
RF_dat <- merge(as.data.frame(random_parms), GSA_summary_dat, by.x = "V36", by.y = "parm_number")
# rename so it's easy to navigate and plot
colnames(RF_dat) <- c("parameter number", "compartments", "event order", 
                      "\u03B2_F,RR", "\u03B2_F,WR", "\u03B2_F,WW", "\u03C3_\u03B2,F", 
                      "\u03B2_D,RR", "\u03B2_D,WR", "\u03B2_D,WW", "\u03C3_\u03B2,D", 
                      "\u03BC_S,RR", "\u03BC_S,WR", "\u03BC_S,WW", "\u03BC_I,RR", "\u03BC_I,WR", "\u03BC_I,WW", "\u03C3_\u03BC,S", "\u03C3_\u03BC,I", 
                      "\u03B3_RR", "\u03B3_WR", "\u03B3_WW", "\u03C3_\u03B3", 
                      "\u03BB_RR", "\u03BB_WR", "\u03BB_WW", "\u03BB_sd", "\u03BD",
                      "K", "\u03C3_K",
                      "t_q", "d", 
                      "Q_0", "I_0", "t_f", "transmission case", colnames(RF_dat[37:55]))
# note the first two variables are categorical
RF_dat$"compartments" <- as.factor(RF_dat$"compartments"); RF_dat$"event order" <- as.factor(RF_dat$"event order"); RF_dat$"transmission case" <- as.factor(RF_dat$"transmission case")

# trying with ratios instead of raw values, e.g., ratio of RR:WW
RF_ratio <- data.frame(compartments = RF_dat$"compartments", 
                       "event order" = RF_dat$"event order", 
                       "allele benefit, enviro. trans. blocking" = RF_dat$"\u03B2_F,RR"/RF_dat$"\u03B2_F,WW", 
                       "wild-type enviro. tranmission" = RF_dat$"\u03B2_F,WW", 
                       "enviro. tranmission variation" = RF_dat$"\u03C3_\u03B2,F",
                       "transmission type" = RF_dat$"transmission case", # want to be able to seperate cases with and without enviro
                       "allele benefit, dens. trans. blocking" = RF_dat$"\u03B2_D,RR"/RF_dat$"\u03B2_D,WW",
                       "wild-type dens. tranmission" = RF_dat$"\u03B2_D,WW",
                       "dens tranmission variation" = RF_dat$"\u03C3_\u03B2,D",
                       "wild-type natural mortality" = RF_dat$"\u03BC_S,WW",
                       "allele benefit, mortality-blocking" = RF_dat$"\u03BC_I,RR"/RF_dat$"\u03BC_I,WW",
                       "wild-type disease mortality" = RF_dat$"\u03BC_I,WW",
                       "natural mortality variation" = RF_dat$"\u03C3_\u03BC,S", 
                       "disease mortality variation" = RF_dat$"\u03C3_\u03BC,I", 
                       "allele benefit, clearance aug." = RF_dat$"\u03B3_RR"/RF_dat$"\u03B3_WW", 
                       "wild-type recovery" = RF_dat$"\u03B3_WW", 
                       "recovery variation" = RF_dat$"\u03C3_\u03B3", 
                       "allele cost, fecundity" = RF_dat$"\u03BB_RR"/RF_dat$"\u03BB_WW",
                       "wild-type fecundity" = RF_dat$"\u03BB_WW"
) # then everything after is the same
RF_ratio <- cbind(RF_ratio, RF_dat[, 27:35]) # , RF_dat[,c(43:44)]) # RF_dat add is for and r/inf info
colnames(RF_ratio) <- c("compartments", "event order",
                        "allele benefit:\nenv. trans. block.", "wild type env. tranmission", "env. tranmission variation", 
                        "transmission type", 
                        "allele benefit:\ndens. trans. block.", "wild type dens. tranmission", "dens. tranmission variation",
                        "wild type natural mortality", "allele benefit:\nmortality block.", "wild type disease mortality", "natural mortality variation", "disease mortality variation",
                        "allele benefit:\nclearance aug.", "wild type recovery", "recovery variation",
                        "allele cost: fecundity", "wild type fecundity", "fecundity variation", 
                        "mutation rate", 
                        "carrying capacity", "carrying capactity variation",
                        "initalization length",
                        "disease gens.\nper host gen.", 
                        "allele frequency after initalization", "initial disease prevelence", "total simulation time")
saveRDS(colnames(RF_ratio), file = "dat/RF_names.Rdata")

# with ratios
forest_extinct_ratio <- randomForest(x=RF_ratio, y=RF_dat$`P(extinct)`, ntree=1000, importance = T, localImp = T)
plot(forest_extinct_ratio) # this is checking the convergence?
varImpPlot(forest_extinct_ratio)
forest_extinct_ratio # checing var explained etc
saveRDS(forest_extinct_ratio, file = "dat/RF_ext_ratio.Rdata")
# ER
forest_ER_ratio <- randomForest(x=RF_ratio, y=RF_dat$`P(ER_50_45)`, ntree=1000, importance = T, localImp = T)
plot(forest_ER_ratio) # this is checking the convergence?
varImpPlot(forest_ER_ratio)
forest_ER_ratio # checking var explained etc
saveRDS(forest_ER_ratio, file = "dat/RF_ER_ratio.Rdata")
# DR -- skip in final fig
forest_DR_ratio <- randomForest(x=RF_ratio, y=RF_dat$`P(DR_50_45)`, ntree=1000, importance = T, localImp = T)
plot(forest_DR_ratio) # this is checking the convergence?
varImpPlot(forest_DR_ratio)
forest_DR_ratio # checking var explained etc
saveRDS(forest_DR_ratio, file = "dat/RF_IL_ratio.Rdata")
# first time to K for ER+
forest_K95_ratio <- randomForest(x=RF_ratio, y=RF_dat$`ER_K95`, ntree=1000, importance = T, localImp = T)
plot(forest_K95_ratio) # this is checking the convergence?
varImpPlot(forest_K95_ratio)
forest_K95_ratio # checking var explained etc -- not super explainatory or surprising....?
saveRDS(forest_K95_ratio, file = "dat/RF_K95_ratio.Rdata")

# plotting non-interacting effect:
partialPlot(forest_ER_ratio, RF_ratio, "event_order")
