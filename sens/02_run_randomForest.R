library(dplyr)
library(randomForest)
library(ggplot2)
library(rpart)
library(rpart.plot)
library(optRF)
library(ranger)
library(cowplot)

# loading the data -- note weird directory because of HPC use
sum_a <- readRDS("dat/gsa_result_0928.Rdata")
sum_b <- readRDS("dat/gsa_result_0930.Rdata")
sum_c <- readRDS("dat/gsa_result_1003.Rdata")
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

# put together
sim_dat <- rbind(a_dtf, b_dtf, c_dtf) # just one???

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

# ranger -- super fast, little visualization, can get importance
frst_clss <- ranger(clss ~ ., data = RF_df, 
                    num.trees = 1000, 
                    classification = T, # make sure it's a classification tree
                    importance = "permutation", # get some importance stats
                    write.forest = F) # save memory
saveRDS(frst_clss, "sens/forest_class.Rdata")

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
p_uk <- RF_class %>% group_by(param_num) %>% summarise(p_uk = sum(clss=="UNK")/1000) # change to n sims
# merge everyone
tmp <- merge(p_ex, p_er, all = T)
tmp <- merge(tmp, p_il, all = T)
tmp <- merge(tmp, p_nr, all = T)
out_mat <- merge(tmp, p_uk, all = T)
# check that out_mat is adding to 1
tmp_check <- out_mat %>% mutate(tot_prob = p_ex + p_er + p_il + p_nr + p_uk); range(tmp_check$tot_prob)
# clean up the input data -- only need one row for each now
pred_vars <- c(1, 2, 3, 6, 10, 14, 17, 22, 26, 28:30, 32, 34, 36, 37, 66, 67, 68, 69)
in_mat <- RF_class[seq(from = 1, to = dim(RF_df)[1], by = 1000), pred_vars]

# now run the forest
tmp_df <- cbind(in_mat[2:20], out_mat[2:6]) 
# consider looking into the quantile option?
clss_qt <- quantreg(cbind(p_ex, p_er, p_il, p_nr, p_uk) ~ ., data = tmp_df, importance = "permute")
# note previous effor to use rfscr... not quite sure the difference
clss_mv <- rfsrc(cbind(p_ex, p_er, p_il, p_nr, p_uk) ~ ., data = tmp_df, importance = "permute")
plot(clss_mv) # first look for one variable only!
imp <- vimp(clss_mv, importance ="permute") # this repeats variable calc from above fyi
# imp$regrOutput$p_ex$importance to get importance for a particular outcome
# look into holdout.vimp -- maybe more what we're looking for? 
ho_imp <- holdout.vimp(cbind(p_ex, p_er, p_il, p_nr, p_uk) ~ ., data = tmp_df)
# also plot.variable -- basically partials for p_er in this case at least
plot.variable(clss_mv, m.target = "p_er") 

# plotting decision tree....
my_tr <- get.tree(clss_mv, tree.id = 1:500, ensemble = T)
plot(my_tr)

# make a data frame of the importances
imp_plt_dt <- data.frame(called_names = imp$xvar.names, 
                         p_ex = imp$regrOutput$p_ex$importance, 
                         p_er = imp$regrOutput$p_er$importance, 
                         p_il = imp$regrOutput$p_il$importance, 
                         p_nr = imp$regrOutput$p_nr$importance, 
                         p_uk = imp$regrOutput$p_uk$importance)
imp_plt_dt$neat_names <- c("compartments", "event order",
                       "enviro. transmission", "dens. transmission",
                       "background mort.", "infect. mort.",
                       "recovery rate",
                       "avg. reproduction",
                       "mutation prob.",
                       "carrying capacity", "carrying capacity SD",
                       "disease gens.",
                       "init. infect.",
                       "transmission type",
                       "adaptation pathway", "adaptive benefit",
                       "dominance", "trait SD", "adaptive cost")
imp_plt_long <- pivot_longer(imp_plt_dt, cols = 2:6)
imp_plt_long$name <- recode(imp_plt_long$name, 
                            p_er = "P(ER)", 
                            p_ex = "P(Ext)", 
                            p_il = "P(IL)", 
                            p_nr = "P(NR)", 
                            p_uk = "Unknown")

inv_imp <- ggplot(data = imp_plt_long, aes(`value`, reorder(neat_names, `value`))) +
  geom_linerange(aes(xmin = 0, xmax = `value`)) +
  geom_point(aes(col = name), size = 3) + # alt: aes(size = log(IncNodePurity))
  labs(x = "", y = NULL, col = "importance to probaiblity:") +
  scale_color_manual(values = c("#ac1457", "black", "#DB6341", "#f1c4a2", "gray")) +
  facet_wrap(~name, nrow = 1) + 
  theme_bw() +
  theme(text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "bottom")

# looking at a tree -- not really helpful???
get.tree(clss_mv, ensemble = T, node.depth = 5)

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
  filter(extinct == 0 & abs(last_50-first_50-tot_50) < 0.5*tot_50 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 1 & final_r_allele > 0.5*max_r_allele ) %>% 
  summarize(`P(ER, Persist.)` = n()/2500)
er_tmp_dat <- tmp %>% 
  group_by(parm_number) %>% 
  filter(extinct == 0 & abs(last_50-first_50-tot_50) < 0.5*tot_50 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 1 & final_inf_prev == 0 & final_r_allele < 0.5*max_r_allele) %>% 
  summarize(`P(ER, Temp.)` = n()/2500)
il_dat <- tmp %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) < 0.5*tot_50 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 0 & final_inf_prev == 0) %>% 
  group_by(parm_number) %>% 
  summarize(`P(IL)` = n()/2500)
# merge
tmp <- merge(prob_dat, er_dat, by = c("parm_number"), all = T)
tmp <- merge(tmp, il_dat, by = c("parm_number"), all = T)
tmp <- merge(tmp, er_tmp_dat, by = c("parm_number"), all = T)
plot_dat <- merge(tmp, cases, by.x = c("parm_number"), by.y = "number")
# rename
plot_dat$compartments <- recode(plot_dat$compartments, "1" = "SIX", "2" = "SIS", "3" = "SIR")
plot_dat$`transmission type` <- recode(plot_dat$`transmission type`, "1" = "density", "2" = "environmental", "3" = "density + environmental")
plot_dat$robustness <- recode(plot_dat$robustness, "1" = "MB", "2" = "TB", "3" = "RA", "4" = "N")
plot_dat$`evolutionary pathway` <- factor(plot_dat$robustness, levels = c("N", "TB", "MB", "RA"))
# convert to long
plot_long <- pivot_longer(plot_dat, cols = 2:5, names_to = "outcome", values_to = "value")
plot_long$value <- ifelse(is.na(plot_long$value), 0, plot_long$value)
plot_long$plus <- plot_long$value + intv
plot_long$less <- plot_long$value - intv
for (i in 1: length(plot_long$value)) {
  if (plot_long$plus[i] > 1) plot_long$plus[i] <- 1
  if (plot_long$less[i] < 0) plot_long$less[i] <- 0
}
plot_long$outcome_f <- factor(plot_long$outcome, levels = c("P(Ext)", "P(ER, Persist.)", "P(ER, Temp.)", "P(IL)"))
# now use predict to predict what those probabilities would be from the RF
base_vect <- c(2, 4, # SIS, BGM -- can comp with 4 for MBG
               -1, -1, -1, 0.05, # 3-6 no enviro dep transmission
               -1, -1, -1, 0.05, # 7-10 always dens dep transmission
               0.05, 0.05, 0.05, # 11-13 + 17 background mort 
               1, 1, 1, # 14-16 + 18 disease mort
               0.05, 0.05, # 17-18 mort sd
               0.05, 0.05, 0.05, 0.05, # 19-22 recovery, unused bc SIX (1)
               1.9, 2, 2.1, 0.05, 0.005, # 23-27 reproduction & mutation
               # 1.7, 1.9, 2.1, 0.05, 0.005, # 23-27 reproduction & mutation
               100, 4, # 28-29 carrying capacity things
               120, 4, # 30-32 timing things # note change from d = 1 to d = 2???
               0.1, 0.1, 160) # init R and init disease, ngens
pred_results_str <- NULL
for (i in 1:N) {
  
  case_vect <- base_vect
  case_vect[35] <- rep_cases$`number`[i]
  # set disease transmission
  if (rep_cases$`transmission type`[i]==1) {
    case_vect[7:9] <- 1.25 # transmission increased (dens)
  }
  if (rep_cases$`transmission type`[i]==2) {
    case_vect[3:5] <- 2 # transmission increased (freq)
  }
  if (rep_cases$robustness[i] == 1) {
    case_vect[14] <- 0
    case_vect[15] <- (case_vect[14]+case_vect[16])/2
  } # 1 = mortality
  if (rep_cases$robustness[i] == 2) {
    if (rep_cases$`transmission type`[i]==1) {case_vect[7:9] <- c(0, 0.625, 1.25)}
    if (rep_cases$`transmission type`[i]==2) {case_vect[3:5] <- c(0, 1, 2)}
  } # 1 = transmission
  if (rep_cases$robustness[i] == 3) {
    case_vect[19] <- 1.516
    case_vect[20] <- (case_vect[19]+case_vect[21])/2
  } # 1 = recovery, nb: recovery higher
  
  # modify base parameters -- compartment
  case_vect[1] <- rep_cases$compartments[i]
  
  pred_results_str <- cbind(pred_results_str, predict(clss_mv, case_vect))
  
}


# points plot
fig1B <- ggplot(plot_long, aes(`evolutionary pathway`, value, col = `transmission type`)) + 
  coord_cartesian(ylim = c(0, NA)) +
  geom_linerange(aes(ymin = less, ymax = plus), lwd = 1) + 
  geom_hline(data = plot_long %>% filter(`evolutionary pathway` == "N"), aes(yintercept = value), col = "gray70", lty = "dashed") + 
  geom_point(size = 2) + scale_color_manual(values = c("#ac1457", "#f1c4a2")) +
  facet_grid(rows = vars(outcome_f), cols = vars(compartments)) + 
  labs(x = "Adaptive pathway", y = "Probability", col = "Transmission type") + 
  theme_bw()  + 
  theme(text = element_text(size = 11), legend.position = "bottom") 