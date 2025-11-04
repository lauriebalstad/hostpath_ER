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

# loading the data -- note weird directory because of HPC use
sum_a <- readRDS("dat/gsa_result_0928.Rdata")
sum_b <- readRDS("dat/gsa_result_0930.Rdata")
sum_c <- readRDS("dat/gsa_result_1003.Rdata")
sum_d <- readRDS("dat/gsa_result_1022.Rdata")
sum_e <- readRDS("dat/gsa_result_1024.Rdata")
sum_f <- readRDS("dat/gsa_result_1025.Rdata")
random_parms <- readRDS("dat/mat_var_0922.Rdata") # figure out what these values were converging around....? new values have super little er
bd_parms <- readRDS("dat/mat_bd_0922.Rdata") # figure out what these values were converging around....? new values have super little er

# convert to df
a_dtf <- as.data.frame(matrix(unlist(sum_a), ncol = 29, byrow = T))
colnames(a_dtf) <- c("extinct", 
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
b_dtf <- as.data.frame(matrix(unlist(sum_b), ncol = 29, byrow = T))
colnames(b_dtf) <- c("extinct",
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

c_dtf <- as.data.frame(matrix(unlist(sum_c), ncol = 29, byrow = T))
colnames(c_dtf) <- c("extinct",
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

d_dtf <- as.data.frame(matrix(unlist(sum_d), ncol = 29, byrow = T))
colnames(d_dtf) <- c("extinct",
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

e_dtf <- as.data.frame(matrix(unlist(sum_e), ncol = 29, byrow = T))
colnames(e_dtf) <- c("extinct",
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

f_dtf <- as.data.frame(matrix(unlist(sum_f), ncol = 29, byrow = T))
colnames(f_dtf) <- c("extinct",
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

# put together
sim_dat <- rbind(a_dtf, b_dtf, c_dtf, d_dtf, e_dtf, f_dtf) # just one???

# change -1s back to 0 bc that's their true rate for random_parms
random_parms <- as.data.frame(random_parms)
for (i in 1:dim(random_parms)[1]) {
  if (all(random_parms[i,3:5] == -1)) random_parms[i,3:5] <- 0 
  if (all(random_parms[i,7:9] == -1)) random_parms[i,7:9] <- 0
}

# trying as a classifier
# set up data: merge parameters and simulated data
RF_class <- merge(random_parms, sim_dat, by.x = "V36", by.y = "parm_number")
RF_class <- merge(RF_class, bd_parms, by.x = "V36", by.y = "number")
# rename things (note need to keep these nice for RF)
colnames(RF_class) <- c("param_num", "compartments", "order", 
                        "B_F_RR", "B_F_WR", "B_F_WW", "SD_B_F", 
                        "B_D_RR", "B_D_WR", "B_D_WW", "SD_B_D", 
                        "M_S_RR", "M_S_WR", "M_S_WW", "M_I_RR", "M_I_WR", "M_I_WW", "SD_M_S", "SD_M_I", 
                        "G_RR", "G_WR", "G_WW", "SD_G", 
                        "L_RR", "L_WR", "L_WW", "SD_L", "mutation",
                        "K", "SD_K",
                        "init_time", "disease_gen", 
                        "init_allele", "init_inif", "final_time", "transmission_cat", "adaptive_cat", colnames(RF_class[38:69]))
# note the the variables that are categorical
RF_class$"compartments" <- as.factor(RF_class$"compartments")
RF_class$"order" <- as.factor(RF_class$"order")
RF_class$"transmission_cat" <- as.factor(RF_class$"transmission_cat")
RF_class$"adaptive_cat" <- as.factor(RF_class$"adaptive_cat")

# figure out classifications
RF_class$ER <- ifelse(RF_class$extinct == 0 & 
                        abs(RF_class$last_50-RF_class$first_50-RF_class$tot_50) < 0.5*RF_class$tot_50 & 
                        RF_class$tot_50 > 3 & 
                        RF_class$at_K95 == 1 & 
                        RF_class$r_allele_peak45 == 1, 
                      1, 0) 

RF_class$IL <- ifelse(RF_class$extinct == 0 & 
                        abs(RF_class$last_50-RF_class$first_50-RF_class$tot_50) < 0.5*RF_class$tot_50 & 
                        RF_class$tot_50 > 3 & 
                        RF_class$at_K95 == 1 & 
                        RF_class$r_allele_peak45 == 0 & 
                        RF_class$final_inf_prev == 0, 
                      1, 0) 

# have extinction column already
RF_class$Ext <- RF_class$extinct 

RF_class$NR <- ifelse(RF_class$extinct == 0 & 
                      (abs(RF_class$last_50-RF_class$first_50-RF_class$tot_50) > 0.5*RF_class$tot_50 | RF_class$tot_50 < 3) &
                      RF_class$at_K == 1, 1, 0) 

RF_class$clss <- ifelse(RF_class$Ext == 1, "Ext", 
                        ifelse(RF_class$ER == 1, "ER", 
                               ifelse(RF_class$IL == 1, "IL",
                                      ifelse(RF_class$NR == 1, "NR", "UNK")))) # ~20% unknowns

RF_class$clss <- as.factor(RF_class$clss) 
# plot(RF_class$clss) # mostly ER, extinction, some NR/UNK

# grab subset (uncorrelated variables) to run RF on
pred_vars <- c(2, 3, 6, 10, 14, 17, 22, 26, 28:30, 32, 34, 36, 37, 66, 67, 68, 69)
RF_df <- cbind(RF_class[, pred_vars], RF_class$clss)
colnames(RF_df) <- c(colnames(RF_df)[1:19], "clss")

# use optRF functions to get a sense of tree optimization
# note this takes a while to run
# exp_opt <- opt_importance(y = RF_df[,20], X=RF_df[,-20])
# saveRDS(exp_opt, "sens/opt_trees.Rdata")
# suggests 1000 trees, tho pretty high stability even at 250

# # ranger -- super fast, little visualization, can get importance
# frst_clss <- ranger(clss ~ ., data = RF_df, 
#                     num.trees = 1000, 
#                     classification = T, # make sure it's a classification tree
#                     importance = "permutation", # get some importance stats
#                     write.forest = F) # save memory
# saveRDS(frst_clss, "sens/forest_class.Rdata")

# randomForest -- slower, but more features??
# rf_clss <- randomForest(clss ~ ., data = RF_df, ntree=50, importance = T, localImp = T, type = "classification")

# rpart for tree
# rp_clss <- rpart(clss ~ ., data = RF_df, method = "class")
# rpart.plot(rp_clss, 
#            type = 5, 
#            legend.x = NA, legend.y = NA,
#            # colors might need reordering, trying to match to decision tree in methods
#            box.palette=list("#ac1457", "black", "#DB6341", "#f1c4a2", "white"),
#            # just play with this... no clear ordering/meaning?
#            col = c("white", "white", "white", "black", "white", "black", "black", "black","black")) # hmmm...
# 
# # plot importance
# RF_imp <- data.frame(name_vals = colnames(RF_df[1:19]), 
#                      typ = c("disease ecology", "disease ecology", 
#                              "force of disease", "force of disease", 
#                              "background", "force of disease", 
#                              "force of disease", 
#                              "background", "adaptation", 
#                              "background", "background", 
#                              "force of disease", 
#                              "background", 
#                              "disease ecology", 
#                              "adaptation", "adaptation", "adaptation", 
#                              "background", 
#                              "adaptation"), 
#                      imprt = RF$variable.importance)
# RF_imp$neat_names <- c("compartments", "event order", 
#                        "enviro. transmission", "dens. transmission", 
#                        "background mort.", "infect. mort.", 
#                        "recovery rate",
#                        "avg. reproduction", 
#                        "mutation prob.", 
#                        "carrying capacity", "carrying capacity SD", 
#                        "disease gens.", 
#                        "init. infect.", 
#                        "transmission type", 
#                        "adaptation pathway", "adaptive benefit", 
#                        "dominance", "trait SD", "adaptive cost")
# 
# RF_imp <- RF_imp %>% arrange(desc(`imprt`))
# RF_plt <- ggplot(data = RF_imp, aes(`imprt`, reorder(neat_names, `imprt`))) +
#   geom_linerange(aes(xmin = 0, xmax = `imprt`)) +
#   geom_point(aes(col = typ), size = 3) + # alt: aes(size = log(IncNodePurity))
#   labs(x = "", y = NULL, col = "parameter\ntype:") +
#   scale_color_manual(values = c("#ac1457", "#DB6341", "#f1c4a2", "black"),
#                      breaks = c("adaptation", "force of disease", "disease ecology", "background"),
#                      labels = c("adaptation", "force\nof disease", "disease\necology", "background")) +
#   theme_bw() +
#   theme(text = element_text(size = 12),
#         legend.text = element_text(size = 12),
#         legend.title = element_text(size = 12),
#         legend.position = "bottom")
# 
# # plot the different outcomes
# RF_dist <- RF_class; RF_dist$clss_lvl <- factor(RF_dist$clss, levels = c("Ext", "ER", "IL", "NR", "UNK"))
# RF_dist_plt <- ggplot(data = RF_dist, aes(clss_lvl, fill = clss_lvl)) + geom_bar(col = "black") + 
#   labs(x = "simulation outcome", y = "count", fill = NULL) + 
#   scale_fill_manual(values = c("black", "#ac1457", "#DB6341", "#f1c4a2", "white")) + 
#   theme_bw() + 
#   theme(text = element_text(size = 12), 
#         legend.position = "none") 

# alternatively, try distributed random forest, which deals with the fact that the outcome is a distribution not a category
# make outcome matrix
# get probabilities
p_ex <- RF_class %>% group_by(param_num) %>% summarise(p_ex = mean(Ext))
p_er <- RF_class %>% group_by(param_num) %>% summarise(p_er = mean(ER))
p_il <- RF_class %>% group_by(param_num) %>% summarise(p_il = mean(IL))
p_nr <- RF_class %>% group_by(param_num) %>% summarise(p_nr = mean(NR))
p_uk <- RF_class %>% group_by(param_num) %>% summarise(p_uk = sum(clss=="UNK")/2500) # change to n sims
# merge everyone
tmp <- merge(p_ex, p_er, all = T)
tmp <- merge(tmp, p_il, all = T)
tmp <- merge(tmp, p_nr, all = T)
out_mat <- merge(tmp, p_uk, all = T)
# check that out_mat is adding to 1
tmp_check <- out_mat %>% mutate(tot_prob = p_ex + p_er + p_il + p_nr + p_uk); range(tmp_check$tot_prob)
# clean up the input data -- only need one row for each now
pred_vars <- c(1, 2, 3, 6, 10, 14, 17, 22, 26, 28:30, 32, 34, 36, 37, 66, 67, 68, 69)
in_mat <- RF_class[seq(from = 1, to = dim(RF_df)[1], by = 2500), pred_vars]

# now run the forest
tmp_df <- cbind(in_mat[2:20], out_mat[2:6]) 
# consider looking into the quantile option?
# clss_qt <- quantreg(cbind(p_ex, p_er, p_il, p_nr, p_uk) ~ ., data = tmp_df, importance = "permute")
# note previous effor to use rfscr... not quite sure the difference
clss_mv <- rfsrc(Multivar(p_ex, p_er, p_il, p_nr, p_uk) ~ . + compartments*adaptive_pathway, data = tmp_df, 
                 importance = "permute", ntree = 500)
plot(clss_mv) # first look for one outcome only!
imp <- vimp(clss_mv, importance ="permute") # this repeats variable calc from above fyi
# use ``imp$regrOutput$p_ex$importance`` to get importance for a particular outcome
# look into holdout.vimp -- maybe more what we're looking for? nope....
# ho_imp <- holdout.vimp(cbind(p_ex, p_er, p_il, p_nr, p_uk) ~ ., data = tmp_df)
# also plot.variable -- basically partials for p_er in this case at least
# plot.variable(clss_mv, m.target = "p_er") 

# plotting decision tree....
my_tr <- get.tree(clss_mv, tree.id = 1:500, ensemble = T)
plot(my_tr)

# looking at a tree -- not really helpful???
# get.tree(clss_mv, ensemble = T, node.depth = 5)

# try predicting a few rows of fig1B, to check out how off base things are
# merge data w/original cases
str_dat <- readRDS("dat/str_fig_1021.Rdata")
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
#   summarize(`P(ER)` = n()/2500)
er_dat <- tmp %>% 
  group_by(parm_number) %>% 
  filter(extinct == 0 & abs(last_50-first_50-tot_50) < 0.5*tot_50 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 1) %>% 
  summarize(`P(ER)` = n()/2500)
# only care about ER as a general category right now
# er_tmp_dat <- tmp %>%  
#   group_by(parm_number) %>% 
#   filter(extinct == 0 & abs(last_50-first_50-tot_50) < 0.5*tot_50 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 1 & final_inf_prev == 0 & final_r_allele < 0.5*max_r_allele) %>% 
#   summarize(`P(ER, Temp.)` = n()/2500)
il_dat <- tmp %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) < 0.5*tot_50 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 0 & final_inf_prev == 0) %>% 
  group_by(parm_number) %>% 
  summarize(`P(IL)` = n()/2500)
# add in a no risk category
nr_dat <- tmp %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) > 0.5*tot_50 & at_K95 == 1) %>% 
  group_by(parm_number) %>% 
  summarize(`P(NR)` = n()/2500)
# merge
tmp <- merge(prob_dat, er_dat, by = c("parm_number"), all = T)
tmp <- merge(tmp, il_dat, by = c("parm_number"), all = T)
tmp <- merge(tmp, nr_dat, by = c("parm_number"), all = T)
plot_dat <- merge(tmp, cases, by.x = c("parm_number"), by.y = "number")

# now use predict to predict what those probabilities would be from the RF
base_vect <- c(2, 4, # 1, 2 -- compartments, event order
               -1, -1, # 3 & 4 -- freq & density transmission rate WW
               0.05, 1, 0.05, 2.1, 0.005, 100, 4, 4, 0.1,
               0, # 14 -- transmission type 
               0, 1, # 15, 16 -- pathway, benefit
               0.5, 0.05, 1.7/2.1)
pred_results_str <- NULL
rep_cases <- cases
# loop through w predict
for (i in 1:N) {
  
  case_vect <- base_vect
  
  # set disease transmission
  if (rep_cases$`transmission type`[i]==1) {
    case_vect[4] <- 1.25 # transmission increased (dens)
    case_vect[3] <- 0 # by def for this figure/RF
    case_vect[14] <- 1
  }
  if (rep_cases$`transmission type`[i]==2) {
    case_vect[3] <- 2 # transmission increased (freq)
    case_vect[4] <- 0
    case_vect[14] <- 2
  }
  if (rep_cases$robustness[i] == 1) {
    case_vect[15] <- 1
    case_vect[16] <- 0 # benefit is defined as the wild type*benefit for robust type (100% reduction here)
  } # 1 = mortality
  if (rep_cases$robustness[i] == 2) {
    case_vect[15] <- 2
    case_vect[16] <- 0 # benefit is defined as the wild type*benefit for robust type (100% reduction here)
  } # 1 = transmission
  if (rep_cases$robustness[i] == 3) {
    case_vect[15] <- 3
    case_vect[16] <- 0.44 # benefit is defined as the wild type*benefit for robust type (100% reduction here)
  } # 1 = recovery, nb: recovery higher
  
  # modify base parameters -- compartment
  case_vect[1] <- rep_cases$compartments[i]
  
  case_vect <- as.data.frame(t(case_vect)); colnames(case_vect) <- clss_mv[["xvar.names"]]
  
  tmp <- predict(clss_mv, case_vect)
  
  out_dat <- c(tmp$regrOutput$p_er$predicted, 
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
colnames(pred_results_str) <- c("parm_number", "P(ER)", "P(Ext)", "P(IL)", "P(NR)", "Unknown")
# double check that things are adding to 1
pred_check <- apply(pred_results_str[, 2:6], 1, sum); range(pred_check)
pred_dat <- merge(pred_results_str, cases, by.x = c("parm_number"), by.y = "number")

# pivot longer for both
pred_dat$shp <- rep("predicted")
plot_dat$shp <- rep("simulated")
plot_long <- pivot_longer(plot_dat, cols = 2:5, names_to = "outcome", values_to = "value")
pred_long <- pivot_longer(pred_dat[, c(1:5, 7:10)], cols = 2:5, names_to = "outcome", values_to = "value") # ignore unks
dat_long <- rbind(plot_long, pred_long)

# rename
dat_long$compartments <- recode(dat_long$compartments, "1" = "SIX", "2" = "SIS", "3" = "SIR")
dat_long$`transmission type` <- recode(dat_long$`transmission type`, "1" = "density", "2" = "environmental", "3" = "density + environmental")
dat_long$robustness <- recode(dat_long$robustness, "1" = "MB", "2" = "TB", "3" = "RA", "4" = "N")
dat_long$`evolutionary pathway` <- factor(dat_long$robustness, levels = c("N", "TB", "MB", "RA"))
dat_long$outcome_f <- factor(dat_long$outcome, levels = c("P(Ext)", "P(ER)", "P(IL)", "P(NR)"))

# points plot -- note RF doesn't know what N is, so drop that from this comp
comp_sim_pred <- ggplot(dat_long %>% filter(`evolutionary pathway` != "N"), aes(`evolutionary pathway`, value, col = `transmission type`, pch = shp)) + 
  coord_cartesian(ylim = c(0, NA)) +
  geom_point(size = 2) + scale_color_manual(values = c("#ac1457", "#f1c4a2")) +
  facet_grid(rows = vars(outcome_f), cols = vars(compartments)) + 
  labs(x = "Adaptive pathway", y = "Probability", col = "Transmission type") + 
  theme_bw()  + 
  theme(text = element_text(size = 11), legend.position = "bottom") 
