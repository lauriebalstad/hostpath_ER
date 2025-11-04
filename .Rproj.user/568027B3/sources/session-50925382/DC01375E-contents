library(gridExtra) # NB: overwrite of combine from randomForest
library(tidyverse)
library(cowplot)
library(ggplot2)
library(viridis)
library(ggh4x)
library(ggparty) 

#-----SUPP: COMP SIMULATIONS-----
# note: run random forest first, that's where sim_dat is from
str_dat_A <- readRDS("dat/gsa_result_1029.Rdata")
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
# parms #1-200

# run 02_run_randomForest.R first, that's the combined data for GSAs
str_dtf_B <- sim_dat %>% filter(parm_number < 201)

# dummy data frame of parameter numbers
parm_nums <- data.frame(parm_number = 1:200)

# calc probs and merge for ext, er, il
str_ex_A <- str_dtf_A %>% filter(extinct == 1) %>% group_by(parm_number) %>% summarise(`P(Ex) - A` = n()/2500)
str_ex_B <- str_dtf_B %>% filter(extinct == 1) %>% group_by(parm_number) %>% summarise(`P(Ex) - B` = n()/2500)
fig_ex <- merge(str_ex_A, str_ex_B, by = "parm_number", all = T) # NB: not keeping parameter row when both are zeros... get those back in!
fig_ex <- merge(fig_ex, parm_nums, all = T)
fig_ex$diff <- ifelse(is.na(fig_ex$`P(Ex) - A`), 0, fig_ex$`P(Ex) - A`) - ifelse(is.na(fig_ex$`P(Ex) - B`), 0, fig_ex$`P(Ex) - B`)

str_er_A <- str_dtf_A %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 &  at_K95 == 1 & r_allele_peak45 == 1) %>% group_by(parm_number) %>% summarise(`P(ER) - A` = n()/2500)
str_er_B <- str_dtf_B %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 &  at_K95 == 1 & r_allele_peak45 == 1) %>% group_by(parm_number) %>% summarise(`P(ER) - B` = n()/2500)
fig_er <- merge(str_er_A, str_er_B, by = "parm_number", all = T) # NB: not keeping parameter row when both are zeros... get those back in!
fig_er <- merge(fig_er, parm_nums, all = T)
fig_er$diff <- ifelse(is.na(fig_er$`P(ER) - A`), 0, fig_er$`P(ER) - A`) - ifelse(is.na(fig_er$`P(ER) - B`), 0, fig_er$`P(ER) - B`)

str_il_A <- str_dtf_A %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 0 & final_inf_prev == 0) %>% group_by(parm_number) %>% summarise(`P(IL) - A` = n()/2500)
str_il_B <- str_dtf_B %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 0 & final_inf_prev == 0) %>% group_by(parm_number) %>% summarise(`P(IL) - B` = n()/2500)
fig_il <- merge(str_il_A, str_il_B, by = "parm_number", all = T) # NB: not keeping parameter row when both are zeros... get those back in!
fig_il <- merge(fig_il, parm_nums, all = T)
fig_il$diff <- ifelse(is.na(fig_il$`P(IL) - A`), 0, fig_il$`P(IL) - A`) - ifelse(is.na(fig_il$`P(IL) - B`), 0, fig_il$`P(IL) - B`)

# combine to plot the things
# add column and pivot longer
# fig_ex_long <- pivot_longer(fig_ex[, c(1, 5:7)], cols = 2:4)
# fig_er_long <- pivot_longer(fig_er[, c(1, 5:7)], cols = 2:4)
# fig_il_long <- pivot_longer(fig_il[, c(1, 5:7)], cols = 2:4)
fig_ex_long <- pivot_longer(fig_ex[, c(1, 4)], cols = 2)
fig_er_long <- pivot_longer(fig_er[, c(1, 4)], cols = 2)
fig_il_long <- pivot_longer(fig_il[, c(1, 4)], cols = 2)
fig_ex_long$metric <- rep("P(Ext)")
fig_er_long$metric <- rep("P(ER)")
fig_il_long$metric <- rep("P(IL)")
fig_long <- rbind(fig_ex_long, fig_er_long, fig_il_long)
fig_long$metric <- factor(fig_long$metric, levels = c("P(Ext)", "P(ER)", "P(IL)"))
# fig_long <- pivot_longer(fig_tot, cols = 2) 

# box plots? -- not as informative to the unimodality
comp_sim_box <- ggplot(fig_long, aes(x = value, fill = metric)) + 
  geom_boxplot(col = "black") + 
  scale_fill_manual(values = c("#ac1457", "#DB6341", "#f1c4a2")) + 
  facet_wrap(~metric, nrow = 3) + 
  labs(x = "Difference across two test simulations") + 
  theme_bw() + theme(legend.position = "none")

# histograms -- this is the move!
comp_sim_hist <- ggplot(fig_long, aes(x = value, fill = metric)) + 
  geom_histogram(col = "black", bins = 30) + 
  scale_fill_manual(values = c("gray13", "#ac1457", "#DB6341")) + # almost color matching the main text, e.g., ex as black
  facet_wrap(~metric, nrow = 3) + 
  labs(x = "Difference across two test simulations", y = NULL) + 
  theme_bw() + theme(legend.position = "none", 
                     axis.text.x = element_text(angle = 45, hjust = 1))

# save figure
# png("figs/figure_plot/comp_sim_hist.png",height=170,width=85,res=400,units='mm')
# print(comp_sim_hist)
# dev.off()
comp_sim_hist

intv <- max(
  quantile(fig_ex_long$value, c(0.025, 0.975)),
  quantile(fig_er_long$value, c(0.025, 0.975)),
  quantile(fig_il_long$value, c(0.025, 0.975))
  )
intv <- round (intv, 2)
# use largest for the intervals, throughotu the main sims figures below

#-----MAIN SIMS: READ IN DATA-----
# str_dat <- readRDS("C:/Users/lauri/AppData/Roaming/MobaXterm/home/evorescue_hostpath_str/dat/str_fig_d2_0203.Rdata")
# str_dat <- readRDS("dat/str_fig_d2_0203.Rdata")
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

# cb_dat <- readRDS("C:/Users/lauri/AppData/Roaming/MobaXterm/home/evorescue_hostpath_str/dat/cb_fig_d2_0203.Rdata")
# cb_dat <- readRDS("dat/cb_fig_d2_0203.Rdata")
cb_dat <- readRDS("dat/cb_fig_1021.Rdata")
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
dc_dat <- readRDS("dat/dc_fig_1022.Rdata")
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
pm_dat <- readRDS("dat/pm_fig_1021.Rdata")
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
  theme(text = element_text(size = 11), legend.position = "bottom") + 
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

# fig1B

# then figure 1B is all the GSA figures--make sure 02_run_randomForest.R is run first!
frst_clss <- readRDS("sens/forest_class.Rdata")
# plot importance
RF_imp <- data.frame(name_vals = colnames(RF_df[1:19]), 
                     typ = c("disease ecology", "disease ecology", 
                             "force of disease", "force of disease", 
                             "background", "force of disease", 
                             "force of disease", 
                             "background", "adaptation", 
                             "background", "background", 
                             "force of disease", 
                             "background", 
                             "disease ecology", 
                             "adaptation", "adaptation", "adaptation", 
                             "background", 
                             "adaptation"), 
                     imprt = frst_clss$variable.importance)
RF_imp$neat_names <- c("compartments", "event order", 
                       "enviro. transmission rate", "dens. transmission rate", 
                       "background mort. rate", "infect. mort. rate", 
                       "recovery rate",
                       "avg. reproduction", 
                       "mutation prob.", 
                       "carrying capacity", "carrying capacity SD", 
                       "disease gens.", 
                       "init. infect.", 
                       "transmission type", 
                       "adaptation pathway", "adaptive benefit", 
                       "dominance", "trait SD", "adaptive cost")

# set up colnames for tree
colnames(RF_df) <- c(RF_imp$neat_names, "clss") # rename so things aren't ugly, rpart deals with this fine compared to ranger

# # now plot importance
# RF_imp <- RF_imp %>% arrange(desc(`imprt`))
# RF_imp$main_name <- rep("Global sensivity importance")
# figS2C <- ggplot(data = RF_imp, aes(`imprt`, reorder(neat_names, `imprt`))) + 
#   geom_linerange(aes(xmin = 0, xmax = `imprt`)) + 
#   geom_point(aes(col = typ), size = 1.5) + # alt: aes(size = log(IncNodePurity))
#   labs(x = NULL, y = NULL, col = "parameter\ntype:") + 
#   scale_color_manual(values = c("#ac1457", "#DB6341", "#f1c4a2", "black"), 
#                      breaks = c("adaptation", "force of disease", "disease ecology", "background"), 
#                      labels = c("adaptation", "force\nof disease", "disease\necology", "background")) + 
#   theme_bw() + facet_wrap(~main_name) + 
#   # guides(color = guide_legend(nrow = 2)) +
#   theme(text = element_text(size = 11), 
#         legend.text = element_text(size = 8),
#         legend.title = element_text(size = 8),
#         legend.position = c(0.75, 0.25)) 
# 
# fig1C <- ggplot(data = RF_imp[1:8,], aes(`imprt`, reorder(neat_names, `imprt`))) + 
#   geom_linerange(aes(xmin = 0, xmax = `imprt`)) + 
#   geom_point(aes(col = typ), size = 3) + # alt: aes(size = log(IncNodePurity))
#   labs(x = NULL, y = NULL, col = "parameter type:") + 
#   scale_color_manual(values = c("#ac1457", "#DB6341", "#f1c4a2", "black"), 
#                      breaks = c("adaptation", "force of disease", "disease ecology", "background"), 
#                      labels = c("adaptation", "force\nof disease", "disease\necology", "background")) + 
#   theme_bw() + facet_wrap(~main_name) + 
#   # guides(color = guide_legend(nrow = 2)) +
#   theme(text = element_text(size = 11), 
#         legend.text = element_text(size = 8),
#         legend.title = element_text(size = 8),
#         legend.position = "bottom") 
# 
# # plot the different outcomes
# RF_dist <- RF_class; RF_dist$clss_lvl <- factor(RF_dist$clss, levels = c("Ext", "ER", "IL", "NR", "UNK"))
# RF_dist$main_name  <- rep("Distribution of simulation outcomes")
# figS2A <- ggplot(data = RF_dist, aes(clss_lvl, fill = clss_lvl)) + geom_bar(col = "black") + 
#   labs(x = "simulation outcome", y = "n simulations from GSA", fill = NULL) + 
#   scale_fill_manual(values = c("black", "#ac1457", "#DB6341", "#f1c4a2", "white")) + 
#   theme_bw() + facet_wrap(~main_name) + 
#   theme(text = element_text(size = 8), 
#         legend.position = "none") 
# 
# # rpart for tree
# rp_clss <- rpart(clss ~ ., data = RF_df, method = "class")
# # rpart.plot(rp_clss, 
# #            type = 5, 
# #            legend.x = NA, legend.y = NA,
# #            # colors might need reordering, trying to match to decision tree in methods
# #            box.palette=list("#ac1457", "black", "#DB6341", "#f1c4a2", "white"),
# #            # just play with this... no clear ordering/meaning?
# #            col = c("white", "white", "white", "black", "white", "black", "black", "black","black")) # hmmm...
# # use ggparty for tree making, will have to add legend and adjust plot theme tho
# rp_clss$splits[, 4] <- round(rp_clss$splits[, 4], 2)
# py_clss <- as.party(rp_clss)
# figS2B <- ggparty(py_clss, horizontal = F, terminal_space = 0.3) +
#   geom_edge() + 
#   # add labels, nudge the label for dens. transmission a bit to the right so it's visable
#   geom_edge_label(size = 2) + 
#   # geom_edge_label(id = c(1:12, 15), size = 2) + geom_edge_label(id = 13:14, size = 2, nudge_x = 0.025) + 
#   geom_node_splitvar(size = 2) +
#   # pass list to gglist containing all ggplot components we want to plot for each
#   geom_node_plot(gglist = list(geom_bar(aes(x = "", fill = clss),
#                                         position = position_fill(), col = "black"),
#                                xlab(NULL), ylab(NULL), 
#                                theme_bw(), theme(text = element_text(size = 6), axis.text.x = element_blank()),
#                                scale_fill_manual(values = c("#ac1457", "black", "#DB6341", "#f1c4a2", "white")))
#   ) + theme(plot.margin=unit(c(0,0,0,0), "mm"))

# make a data frame of the importances
imp_plt_dt <- data.frame(called_names = imp$xvar.names,
                         p_ex = imp$regrOutput$p_ex$importance,
                         p_er = imp$regrOutput$p_er$importance,
                         p_il = imp$regrOutput$p_il$importance,
                         p_nr = imp$regrOutput$p_nr$importance,
                         p_uk = imp$regrOutput$p_uk$importance)
imp_plt_dt$neat_names <- c("host response potential", "event order",
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

imp_plt_short <- imp_plt_long %>% 
  filter(neat_names %in% c("host response potential", 
                           "event order",
                           "disease gens.",
                           "transmission type",
                           "adaptation pathway", 
                           "dominance"))

fig1C <- ggplot(data = imp_plt_short, aes(`value`, reorder(neat_names, `value`))) +
  geom_linerange(aes(xmin = 0, xmax = `value`)) +
  geom_point(aes(col = name), size = 3) + # alt: aes(size = log(IncNodePurity))
  labs(x = NULL, y = NULL, col = "importance to event:") +
  scale_color_manual(values = c("#ac1457", "black", "#DB6341", "#f1c4a2", "gray")) +
  facet_wrap(~name, nrow = 1) +
  theme_bw() +
  theme(text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "bottom")

fig1_h <- plot_grid(fig1A, fig1B, ncol = 2, rel_widths = c(0.6, 1))
fig1 <- plot_grid(fig1_h, fig1C, ncol = 1, rel_heights = c(1, 0.5))
# fig1

# combine figure 1
# fig1 <- plot_grid(fig1A, fig1B, rel_widths= c(0.55, 1), nrow = 1)

# save figure
png("figs/figure_plot/pathways_trees_1022.png",height=190,width=170,res=400,units='mm')
print(fig1)
dev.off()

# figure out pieces for supp figure???
# combine the pieces for the supplementary figure
# figS2_i <- ggdraw() +
#   draw_plot(figS2C) +
#   draw_plot(figS2A, x = 0.6, y = 0.3, width = .35, height = .35)
# con_tab <- data.frame("outcome" = c("ER", "Ext", "IL", "NR", "UNK"), 
#                  ER = frst_clss$confusion.matrix[, 1],
#                  Ext = frst_clss$confusion.matrix[, 2],
#                  IL = frst_clss$confusion.matrix[, 3],
#                  NR = frst_clss$confusion.matrix[, 4],
#                  UNK = frst_clss$confusion.matrix[, 5]
# )
# con_grob <- tableGrob(con_tab, rows = NULL, theme = ttheme_default(base_size = 6, parse = T))

# figS2_v <- plot_grid(figS2A, con_grob, ncol = 1, rel_heights = c(1, 1))
# figS2_h <- plot_grid(figS2_v, figS2C, ncol = 2, rel_widths = c(0.7, 1))
# figS2 <- plot_grid(figS2_h, figS2B, ncol = 1, rel_heights = c(1, 1.2))
# # figS2

# png("figs/figure_plot/SUPP_decision_tree_1022.png",height=210,width=170,res=400,units='mm')
# print(figS2)
# dev.off()

#-----FIG2: COST/BENEFIT-----
# merge data: use the percents directly to avoid flipping around later
sis_cb_A <- expand.grid("benefit" = c(0, 50, 75, 100), # writing in prop benefit
                        "cost" = c(1.7, 1.9, 2.1), # 1.7 = some cost, 2.1 = no cost
                        "case" = c("Transmission-blocking"),
                        "compartments" = c(2))
# "compartments" = c(2))
sis_cb_B <- expand.grid("benefit" = c(0, 50, 200, 1000), 
                        "cost" = c(1.7, 1.9, 2.1), # 1.7 = some cost, 2.1 = no cost
                        "case" = c("Recovery-augmenting"),
                        "compartments" = c(2))
sis_cb_C <- expand.grid("benefit" = c(0, 50, 75, 100), 
                        "cost" = c(1.7, 1.9, 2.1), # 1.7 = some cost, 2.1 = no cost
                        "case" = c("Mortality-blocking"),
                        "compartments" = c(2))
sis_cb <- rbind(sis_cb_A, sis_cb_B, sis_cb_C)
N <- dim(sis_cb)[1]
sis_cb$number <- 1:N
cost_ben <- merge(cb_dtf, sis_cb, by.x = "parm_number", by.y = "number", all = T) # %>% filter(cost > 1.8)
# get summaries to plot -- P(ER) + T_K + max/final R/inf (6) for the ER cases only
ex_dat <- cost_ben %>% filter(extinct == 1) %>% 
  group_by(benefit, cost, case, compartments) %>% 
  summarize(`P(Ext)` = n()/2500)
er_dat <- cost_ben %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) < 0.5*tot_50 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 1) %>% 
  group_by(benefit, cost, case, compartments) %>% 
  summarize(`P(ER)` = n()/2500)
dr_dat <- cost_ben %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) < 0.5*tot_50 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 0 & final_inf_prev == 0) %>%    
  group_by(benefit, cost, case, compartments) %>% 
  summarize(`P(IL)` = n()/2500)
plt_dat <- merge(ex_dat, er_dat)
plt_dat <- merge(plt_dat, dr_dat)
# rbind creates a ton of NAs but it's okay...
cb_summ_long <- pivot_longer(plt_dat, cols= 5:7, names_to = "outcome", values_to = "value")
cb_summ_long$plus <- cb_summ_long$value + intv
cb_summ_long$less <- cb_summ_long$value - intv
for (i in 1:length(cb_summ_long$value)) {
  if (cb_summ_long$plus[i] > 1) cb_summ_long$plus[i] <- 1
  if (cb_summ_long$less[i] < 0) cb_summ_long$less[i] <- 0
}
cost_ben$firstK95 <- (cost_ben$firstK95-30)/4 # remove 100 year burn in period (100/2 disease gens)
cost_ben <- rename(cost_ben, 
                 `Final M allele\nfrequency (ER)` = final_r_allele, 
                 `Final infection\nprevelence (ER)` = final_inf_prev, 
                 `Time to\nrecovery (ER)` = firstK95)
cb_all_long <- pivot_longer(cost_ben %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) < 0.5*tot_50 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 1) %>% 
                              select(c(28, 30:33)), 
                            cols = 1, names_to = "outcome", values_to = "value")
cb_summ_long$outcome_f <- factor(cb_summ_long$outcome, levels = c("P(Ext)", "P(ER)", "P(IL)", 
                                                            "Time to\nrecovery (ER)"))
cb_all_long$outcome_f <- factor(cb_all_long$outcome, levels = c("P(Ext)", "P(ER)", "P(IL)", 
                                                            "Time to\nrecovery (ER)"))
# cb_summ_long$percent_benefit_f <- factor(cb_summ_long$percent_benefit, levels = c(0, 0.5, 0.75, 1))
# cb_all_long$percent_benefit_f <- factor(cb_all_long$percent_benefit, levels = c(0, 0.5, 0.75, 1))
# remove the no-evo cases from box plots
# cb_mu_beta <- cb_all_long %>% filter(case %in% c("Mortality-blocking", "Transmission-blocking")) %>% filter(benefit < 1)
# cb_gamma <- cb_all_long %>% filter(case == "Recovery-augmenting") %>% filter(benefit > 0)
# cb_no_evo <- rbind(cb_mu_beta, cb_gamma)
cb_no_evo <- cb_all_long %>% filter(benefit > 0)
cb_no_evo_HL <- cb_no_evo # %>% filter(cost %in% c(1.7, 2.1))
cb_summ_long_HL <- cb_summ_long # %>% filter(cost %in% c(1.7, 2.1))
fig2 <- ggplot(NULL, aes(x=as.factor(benefit), y=value, col = as.factor(cost))) + 
  coord_cartesian(ylim = c(0, NA)) +
  geom_point(data = cb_summ_long_HL, size = 2) + 
  geom_boxplot(data = cb_no_evo_HL) + 
  geom_linerange(data = cb_summ_long_HL, aes(ymin = less, ymax = plus), lwd = 1) + 
  scale_color_manual(values = c("#ac1457", "#DB6341", "#f1c4a2")) +
  # scale_color_manual(values = c("#ac1457", "#f1c4a2")) + 
  facet_grid(rows = vars(outcome_f), cols = vars(case), scales = "free") + 
  theme_bw() + theme(text = element_text(size = 11), legend.position = "bottom") + 
  labs(x = "Mutant allele benefit\n(percent change in probability)", col = "Mutant allele\ncost (avg. fecundity)", y = NULL)
# save figure
png("figs/figure_plot/cost_benefit_1022.png",height=160,width=170,res=400,units='mm')
print(fig2)
dev.off()
# fig2

#-----FIG3: POP SIZES-----
# merge data
sis_pop_mut <- expand.grid("population size" = c(50, 100, 500), 
                           "robustness" = c(1), # M, G only bc they are strongest ER from past exploration
                           "mutation rate" = c(0.0005, 0.005, 0.01),
                           "compartments" = c(2), 
                           "transmission type" = c(1, 2)) 
N <- dim(sis_pop_mut)[1]
sis_pop_mut$number <- 1:N
ex_risk <- merge(pm_dtf, sis_pop_mut, by.x = "parm_number", by.y = "number", all = T)
ex_risk <- ex_risk %>% filter(compartments == 2, robustness == 1, `transmission type`==1)
# rename
ex_risk$robustness <- recode(ex_risk$robustness, "1" = "Mortality-blocking", "2" = "Transmission-blocking", "3" = "Recovery-augementing", "4" = "N")
ex_risk$compartments <- recode(ex_risk$compartments, "1" = "SIX", "2" = "SIS", "3" = "SIR")
ex_risk$transmission <- recode(ex_risk$`transmission type`, "1" = "Density-dependent", "2" = "Environmental")
# get summaries to plot -- P(ER) + T_K + max/final R/inf (6) for the ER cases only
ex_dat <- ex_risk %>% filter(extinct == 1) %>% 
  group_by(parm_number) %>% 
  summarize(`P(Ext)` = n()/2500)
er_dat <- ex_risk %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) < 0.5*tot_50 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 1) %>% 
  group_by(parm_number) %>% 
  summarize(`P(ER)` = n()/2500)
dr_dat <- ex_risk %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) < 0.5*tot_50 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 0 & final_inf_prev == 0) %>%    
  group_by(parm_number) %>% 
  summarize(`P(IL)` = n()/2500)
plt_dat <- merge(ex_dat, er_dat)
plt_dat <- merge(plt_dat, dr_dat)
tmp <- merge(plt_dat, sis_pop_mut, by.x = "parm_number", by.y = "number")
tmp$compartments <- recode(tmp$compartments, "1" = "SIX", "2" = "SIS", "3" = "SIR")
tmp$robustness <- recode(tmp$robustness, "1" = "Mortality-blocking", "2" = "Transmission-blocking", "3" = "Recovery-augementing", "4" = "N")
tmp$transmission <- recode(tmp$`transmission type`, "1" = "Density-dependent", "2" = "Environmental")
# remove NAs --> those are 0s
tmp[is.na(tmp)] <- 0
pm_summ_long <- pivot_longer(tmp[, c(2:8, 10)], cols= 1:3, names_to = "outcome", values_to = "value")
pm_summ_long$plus <- pm_summ_long$value + intv
pm_summ_long$less <- pm_summ_long$value - intv
for (i in 1:length(pm_summ_long$value)) {
  if (pm_summ_long$plus[i] > 1) pm_summ_long$plus[i] <- 1
  if (pm_summ_long$less[i] < 0) pm_summ_long$less[i] <- 0
}
ex_risk$firstK95 <- (ex_risk$firstK95-30)/4 # remove 50 year burn in period
ex_risk <- rename(ex_risk, 
                   `Time to\nrecovery (ER)` = firstK95)
pm_all_long <- pivot_longer(ex_risk %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) < 0.5*tot_50 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 1) %>% 
                              select(c(28, 30:33, 35)), 
                            cols = 1, names_to = "outcome", values_to = "value")
custom_y <- list(
  scale_y_continuous(limits = c(0, 1)),
  scale_y_continuous(limits = c(0, 1)),
  scale_y_continuous(limits = c(0, 140))
)
pm_summ_long$outcome_f <- factor(pm_summ_long$outcome, levels = c("P(Ext)", "P(ER)", "P(IL)", 
                                                                  "Final M allele\nfrequency (ER)", "Final infection\nprevelence (ER)", "Time to\nrecovery (ER)"))
pm_all_long$outcome_f <- factor(pm_all_long$outcome, levels = c("P(Ext)", "P(ER)", "P(IL)", 
                                                                "Final M allele\nfrequency (ER)", "Final infection\nprevelence (ER)", "Time to\nrecovery (ER)"))
fig3 <- ggplot(NULL, aes(x=as.factor(`population size`), y=value, col=as.factor(`mutation rate`))) + 
  geom_boxplot(data = pm_all_long) +
  coord_cartesian(ylim = c(0, NA)) +
  geom_linerange(data = pm_summ_long, aes(ymin = less, ymax = plus), lwd = 1) + 
  geom_point(data = pm_summ_long, size = 2) + 
  scale_color_manual(values = c("#ac1457","#DB6341", "#f1c4a2")) + 
  facet_grid(outcome_f ~ robustness, scales = "free_y") + 
  theme_bw() + theme(text = element_text(size = 11), legend.position = "bottom") + 
  labs(x = "Population size", col = "Mutation rate", y = NULL)
# save figure
png("figs/figure_plot/pop_mut_1022.png",height=160,width=85,res=400,units='mm')
print(fig3)
dev.off()
# if 3x3 grid, width = 170 (this is pm_fig_A_0209 or pm_fig_B_0209, using mortality blocking (2))
# fig3

#-----SUPP: DISEASE CYCLES-----
sis_dc <- expand.grid("disease cycles" = 2:5, # probably provides enough insight...
                      "robustness" = c(1, 3), # M, G only bc they are strongest ER (from fig 3)
                      "event order" = c(2, 4, 6),
                      "transmission" = c(1),
                      "compartments" = c(2)) # 1 = density, 2 = envrionemntal
N <- dim(sis_dc)[1]
sis_dc$parm_number <- 1:N
dc_dat <- merge(dc_dtf, sis_dc, by = "parm_number", all = T)
dc_dat <- dc_dat %>% filter(transmission == 1 & compartments == 2)
sis_dc <- sis_dc %>% filter(transmission == 1 & compartments == 2)
# dc_dat <- dc_dat %>% filter(robustness == 3 & compartments == 2)
# sis_dc <- sis_dc %>% filter(robustness == 3 & compartments == 2)
# rename
dc_dat$`event order` <- recode(dc_dat$`event order`, "2" = "Transmission, Recovery,\nMortality", "4" = "Mortality, Transmission,\nRecovery", "6" = "Recovery, Mortality,\nTransmission")
dc_dat$robustness <- recode(dc_dat$robustness, "1" = "MB", "2" = "TB", "3" = "RA", "4" = "N")
dc_dat$compartments <- recode(dc_dat$compartments, "1" = "SIX", "2" = "SIS", "3" = "SIR")
dc_dat$transmission <- recode(dc_dat$transmission, "1" = "Density-dependent", "2" = "Environmental")
sis_dc$`event order` <- recode(sis_dc$`event order`, "2" = "Transmission, Recovery,\nMortality", "4" = "Mortality, Transmission,\nRecovery", "6" = "Recovery, Mortality,\nTransmission")
sis_dc$robustness <- recode(sis_dc$robustness, "1" = "MB", "2" = "TB", "3" = "RA", "4" = "N")
sis_dc$compartments <- recode(sis_dc$compartments, "1" = "SIX", "2" = "SIS", "3" = "SIR")
sis_dc$transmission <- recode(sis_dc$transmission, "1" = "Density-dependent", "2" = "Environmental")
# get summaries to plot -- P(ER) + T_K + max/final R/inf (6) for the ER cases only --> note more generous def of pop drop
er_dat <- dc_dat %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) < 0.5*tot_50 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 1) %>% 
  group_by(`parm_number`) %>% 
  summarize(`P(ER)` = n()/2500)
tmp <- merge(er_dat, sis_dc, by = "parm_number", all = T)
er_dat <- tmp; er_dat[is.na(er_dat)] <- 0
dr_dat <- dc_dat %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) < 0.5*tot_50 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 0 & final_inf_prev == 0) %>%    
  group_by(parm_number) %>% 
  summarize(`P(IL)` = n()/2500)
tmp <- merge(dr_dat, sis_dc, by = "parm_number", all = T)
dr_dat <- tmp; dr_dat[is.na(dr_dat)] <- 0
ex_dat <- dc_dat %>% filter(extinct == 1) %>% 
  group_by(parm_number) %>% 
  summarize(`P(Ext)` = n()/2500)
tmp <- merge(ex_dat, sis_dc, by = "parm_number", all = T)
ex_dat <- tmp; ex_dat[is.na(ex_dat)] <- 0
# time to K
dc_dat$firstK95 <- (dc_dat$firstK95-30)/4
dc_dat <- rename(dc_dat, 
                   `Final M allele\nfrequency (ER)` = final_r_allele, 
                   `Final infection\nprevelence (ER)` = final_inf_prev, 
                   `Time to\nrecovery (ER)` = firstK95)
dc_timeK_long <- pivot_longer(dc_dat %>% 
                                filter(extinct == 0 & abs(last_50-first_50-tot_50) < 0.5*tot_50 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 1) %>% 
                                select(c(9, 11, 28, 30:34)), 
                              cols = 1:3, names_to = "outcome", values_to = "value")
tmp_dat <- merge(er_dat, ex_dat)
tmp_dat <- merge(tmp_dat, dr_dat)
dc_all_long <- pivot_longer(tmp_dat[, 2:9], cols= 6:8, names_to = "outcome", values_to = "value")
dc_all_long$plus <- dc_all_long$value + intv
dc_all_long$less <- dc_all_long$value - intv
for (i in 1:length(dc_all_long$value)) {
  if (dc_all_long$plus[i] > 1) dc_all_long$plus[i] <- 1
  if (dc_all_long$less[i] < 0) dc_all_long$less[i] <- 0
}
dc_timeK_long$outcome_f <- factor(dc_timeK_long$outcome, levels = c("P(Ext)", "P(ER)", "P(IL)", 
                                                                  "Final M allele\nfrequency (ER)", "Final infection\nprevelence (ER)", "Time to\nrecovery (ER)"))
dc_all_long$outcome_f <- factor(dc_all_long$outcome, levels = c("P(Ext)", "P(ER)", "P(IL)", 
                                                                "Final M allele\nfrequency (ER)", "Final infection\nprevelence (ER)", "Time to\nrecovery (ER)"))
figS1 <- ggplot(NULL, aes(x=as.factor(`disease cycles`), y=value, 
                                  col=as.factor(robustness))) + 
  coord_cartesian(ylim = c(0, NA)) +
  geom_boxplot(data = dc_timeK_long) + 
  geom_linerange(data = dc_all_long, aes(ymin = less, ymax = plus), lwd = 1) + 
  geom_point(data = dc_all_long, size = 2) + 
  scale_color_manual(values = c("#ac1457", "#f1c4a2")) + 
  facet_grid(row = vars(outcome_f), col = vars(`event order`), scales = "free_y") + 
  theme_bw() + theme(text = element_text(size = 11), legend.position = "bottom") + 
  labs(x = "Disease cycles between reproduction", col = "Adaptive pathway", y = NULL)
# check labels/facet options! can either filter by one of the robustness or by one of the transmission 
# save figure
png("figs/figure_plot/dc_supp_1022.png",height=210,width=170,res=400,units='mm')
print(figS1)
dev.off()
# figS1



