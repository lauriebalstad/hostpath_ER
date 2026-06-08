library(gridExtra) 
library(tidyverse)
library(cowplot)
library(ggplot2)
library(viridis)
library(ggh4x)
# library(ggparty) 
# source("sens/02_run_randomForest.R") # NB: only need through RF_class creation

#-----SUPP: COMP SIMULATIONS-----
str_dat_A <- readRDS("dat/gsa_result_0220-1006.Rdata")
# convert to df
str_dtf_A <- as.data.frame(matrix(unlist(str_dat_A), ncol = 22, byrow = T))
colnames(str_dtf_A) <- c("extinct", 
                         "pop_drop20", "pop_drop50", # "pop_drop80",
                         "r_allele_peak15", "r_allele_peak45", #"r_allele_peak75",
                         "final_r_allele", 
                         "final_pop_size", "final_inf_prev",
                         "max_r_allele", # "time_max_r_allele",
                         "max_inf_prev", "time_last_zero_inf",
                         "min_pop", "time_min_pop",
                         "first_20", "first_50", # "first_80",
                         "last_20", "last_50", # "last_80",
                         "tot_20", "tot_50", #"to t_80",
                         "at_K95", "firstK95",
                         # "r_ts_d0", 
                         "parm_number")

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

str_il_A <- str_dtf_A %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 0 & final_inf_prev == 0) %>% group_by(parm_number) %>% summarise(`P(Inf. Loss) - A` = n()/2500)
str_il_B <- str_dtf_B %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) <3 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 0 & final_inf_prev == 0) %>% group_by(parm_number) %>% summarise(`P(Inf. Loss) - B` = n()/2500)
fig_il <- merge(str_il_A, str_il_B, by = "parm_number", all = T) # NB: not keeping parameter row when both are zeros... get those back in!
fig_il <- merge(fig_il, parm_nums, all = T)
fig_il$diff <- ifelse(is.na(fig_il$`P(Inf. Loss) - A`), 0, fig_il$`P(Inf. Loss) - A`) - ifelse(is.na(fig_il$`P(Inf. Loss) - B`), 0, fig_il$`P(Inf. Loss) - B`)

# combine to plot the things
# add column and pivot longer
fig_ex_long <- pivot_longer(fig_ex[, c(1, 4)], cols = 2)
fig_er_long <- pivot_longer(fig_er[, c(1, 4)], cols = 2)
fig_il_long <- pivot_longer(fig_il[, c(1, 4)], cols = 2)
fig_ex_long$metric <- rep("P(Ext)")
fig_er_long$metric <- rep("P(All ER)")
fig_il_long$metric <- rep("P(Inf. Loss)")
fig_long <- rbind(fig_ex_long, fig_er_long, fig_il_long)
fig_long$metric <- factor(fig_long$metric, levels = c("P(Ext)", "P(All ER)", "P(Inf. Loss)"))
# fig_long <- pivot_longer(fig_tot, cols = 2) 

# histograms -- this is the move!
comp_sim_hist <- ggplot(fig_long, aes(x = value, fill = metric)) + 
  geom_histogram(col = "black", bins = 30) + 
  scale_fill_manual(values = c("gray13", "#ac1457", "#DB6341")) + # almost color matching the main text, e.g., ex as black
  facet_wrap(~metric, nrow = 3) + 
  labs(x = "Difference across two test simulations", y = NULL) + 
  theme_bw() + theme(legend.position = "none", 
                     axis.text.x = element_text(angle = 45, hjust = 1))

# save figure
png("figs/figure_plot/SUPP_conv_0220.png",height=170,width=85,res=400,units='mm')
print(comp_sim_hist)
dev.off()

# comp_sim_hist

intv <- max(
  quantile(fig_ex_long$value, c(0.025, 0.975)),
  quantile(fig_er_long$value, c(0.025, 0.975)),
  quantile(fig_il_long$value, c(0.025, 0.975))
  )
intv <- round(intv, 3)
# use largest for the intervals, throughotu the main sims figures below

#-----MAIN SIMS: READ IN DATA-----
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

# cb_dat <- readRDS("dat/cb_fig_0128.Rdata")
cb_dat <- readRDS("dat/cb_fig_0606.Rdata")
# cb_dtf <- as.data.frame(matrix(unlist(cb_dat), ncol = 22, byrow = T))
cb_dtf <- as.data.frame(matrix(unlist(cb_dat), ncol = 27, byrow = T))
colnames(cb_dtf) <- c("extinct", 
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

# dc_dat <- readRDS("dat/dc_fig_0128.Rdata")
dc_dat <- readRDS("dat/dc_fig_0605.Rdata")
dc_dtf <- as.data.frame(matrix(unlist(dc_dat), ncol = 27, byrow = T))
colnames(dc_dtf) <- c("extinct", 
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

# pm_dat <- readRDS("C:/Users/lauri/AppData/Roaming/MobaXterm/home/evorescue_hostpath_str/dat/pm_fig_d2_0203.Rdata")
# pm_dat <- readRDS("dat/pm_fig_d2_0203.Rdata")
pm_dat <- readRDS("dat/pm_fig_0605.Rdata")
pm_dtf <- as.data.frame(matrix(unlist(pm_dat), ncol = 27, byrow = T))
colnames(pm_dtf) <- c("extinct", 
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

#-----FIG1: TIME SERIES & STR-----
# time series data loading
# str_dat <- readRDS("dat/all_ts_dat_0120.rds") # conintue to play w... there's a lot of cases w/o pop drop getting hit
str_dat <- readRDS("dat/all_ts_dat_0605.rds") # conintue to play w... there's a lot of cases w/o pop drop getting hit
# str_dat <- all_ts_dat
# rotate
all_ts_long <- pivot_longer(str_dat %>% select(c(1:2, 8, 9, 11)), cols = c(3:5), names_to = "outcome", values_to = "value")
all_ts_long$outcome_f <- factor(all_ts_long$outcome, 
                                levels = c("Scaled\npopulation size", "Adaptive allele\nfrequency", "Infection\nprevelence"), 
                                labels = c("Scaled pop. size", "Adaptive allele freq.", "Infection prev."))
colors <- c("Extirpation\n(Ext)" = alpha("#ac1457", 0.8), "Persist. ER" = alpha("#DB6341",0.8),
            "Temp. ER" = "#D71B33", "Inf. Loss (IL)" = alpha("#f1c4a2", 0.8),
            "Lasting\ndisease\nshock" = "black")
fig1A <- ggplot(all_ts_long %>% filter(ts > -8, ts < 35), aes(ts, value)) + # 7, 8, 9, 10 gives exinction
  geom_line(aes(group = trial), col = alpha("black", 0.12), lwd = 0.75/1.5) +
  geom_line(data = all_ts_long %>% filter(trial == 12, ts > -8, ts < 35), aes(color = "Lasting\ndisease\nshock"), lwd = 0.8/1.5) + # extinction
  geom_line(data = all_ts_long %>% filter(trial == 3, ts > -8, ts < 35), aes(color = "Persist. ER"), lwd = 0.9/1.5) + # ER
  geom_line(data = all_ts_long %>% filter(trial == 7, ts > -8, ts < 35), aes(color = "Temp. ER"), lwd = 0.9/1.5) + # back selection
  geom_line(data = all_ts_long %>% filter(trial == 10, ts > -8, ts < 35), aes(color = "Inf. Loss (IL)"), lwd = 0.95/1.5) + # DR
  geom_line(data = all_ts_long %>% filter(trial == 5, ts > -8, ts < 35), aes(color = "Extirpation\n(Ext)"), lwd = 0.9/1.5) + # extinction
  facet_wrap(~outcome_f, scale = "free_y", ncol = 1) + theme_bw() + 
  scale_color_manual(values = colors, 
                     breaks = c("Lasting\ndisease\nshock", "Temp. ER", "Inf. Loss (IL)", 
                                "Extirpation\n(Ext)", "Persist. ER")) + 
  theme_bw(base_size = 9.5) +
  theme(legend.position = "bottom", 
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.key.spacing.y = unit(0, "pt")) + 
  guides(color = guide_legend(nrow = 3)) +
  labs(x = "Standardized generation", y = NULL, col = NULL)

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
  summarize(`P(All ER)` = n()/2500, 
            `Time to\nrecovery (ER)` = mean((firstK95-120/3)/3), 
            upp_qnt = quantile((firstK95-120/3)/3, probs = 0.975),
            lwr_qnt = quantile((firstK95-120/3)/3, probs = 0.025))
# # note temporary ER also possible--combine for simplicity
# er_tmp_dat <- str_dtf %>% 
#   group_by(parm_number) %>%
#   filter(extinct == 0 & tot_50 > 0.75*abs(last_50-first_50 + 1) & at_K95 == 1 & r_allele_peak45 == 1  & final_r_allele < 0.25 & final_inf_prev < 0.1) %>%
#   summarize(`P(T-ER)` = n()/2500)
il_dat <- str_dtf %>% filter(extinct == 0 & abs(last_50-first_50-tot_50+1) < 0.8*tot_50 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 0 & final_inf_prev == 0) %>% 
  group_by(parm_number) %>% 
  summarize(`P(IL)` = n()/2500)
# merge summaries together
tmp <- merge(prob_dat, er_dat, by = c("parm_number"), all = T)
tmp <- merge(tmp, il_dat, by = c("parm_number"), all = T)
# tmp <- merge(tmp, er_tmp_dat, by = c("parm_number"), all = T)

# re-combine w cases data
plot_dat <- merge(tmp, cases, by.x = c("parm_number"), by.y = "number")
plot_dat$`Time to\nrecovery (ER)` <- ifelse(plot_dat$robustness == "N", NA, plot_dat$`Time to\nrecovery (ER)`)
plot_dat$upp_qnt <- ifelse(plot_dat$robustness == "N", NA, plot_dat$upp_qnt)
plot_dat$lwr_qnt <- ifelse(plot_dat$robustness == "N", NA, plot_dat$lwr_qnt)
# reorder
plot_dat$compartments <- factor(plot_dat$compartments, levels = c("Mortality", "Recovery", "Immunity"))
plot_dat$`evolutionary pathway` <- factor(plot_dat$robustness, levels = c("N", "TB", "MB", "RA"))
# convert to long
plot_long <- pivot_longer(plot_dat, cols = c(2:4, 7), names_to = "outcome", values_to = "value")
# add the line range information
plot_long$value <- ifelse(is.na(plot_long$value) & plot_long$outcome != "Time to\nrecovery (ER)", 0, plot_long$value)
plot_long$plus <- ifelse(plot_long$outcome == "Time to\nrecovery (ER)", plot_long$upp_qnt, plot_long$value + intv)
plot_long$less <- ifelse(plot_long$outcome == "Time to\nrecovery (ER)", plot_long$lwr_qnt, plot_long$value - intv)
# making sure it's confined to probability space when relevent
for (i in 1:length(plot_long$value)) {
  if (plot_long$plus[i] > 1 & plot_long$outcome[i] != "Time to\nrecovery (ER)") plot_long$plus[i] <- 1
  if (plot_long$less[i] < 0 & plot_long$outcome[i] != "Time to\nrecovery (ER)") plot_long$less[i] <- 0
}
plot_long$outcome_f <- factor(plot_long$outcome, levels = c("P(Ext)", "P(All ER)", "P(IL)", "Time to\nrecovery (ER)"))

# points plot
fig1B <- ggplot(plot_long, aes(`evolutionary pathway`, col = outcome_f, shape = `transmission`)) + 
  geom_linerange(aes(ymin = less, ymax = plus), lwd = 0.65, position = position_dodge(width = 0.7)) + 
  geom_hline(data = plot_long %>% filter(`evolutionary pathway` == "N"), aes(yintercept = value), col = "gray70", lty = "dashed") + 
  geom_point(aes(y = value), size = 1.4, fill = "white", position = position_dodge(width = 0.7)) + 
  scale_shape_manual(values = c(15, 16, 21)) + 
  scale_color_manual(values = c("#ac1457", "#DB6341", "#f1c4a2", "black")) +
  facet_grid(rows = vars(outcome_f), cols = vars(compartments), scales = "free") + 
  labs(x = "Adaptive pathway", y = NULL, shape = NULL) + 
  theme_bw(base_size = 9.5) + guides(col = "none") + 
  theme(legend.position = "bottom", 
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.key.spacing.y = unit(0, "pt"))  
  
# fig1B

# then figure 1C is all the GSA figures--make sure 02_run_randomForest.R is run first!
frst <- readRDS("sens/rf_output.Rdata")

# make a data frame of the importances
imp_plt_dt <- data.frame(called_names = frst$xvar.names,
                         p_ex = frst$regrOutput$p_ex$importance,
                         p_er = frst$regrOutput$p_er$importance,
                         p_il = frst$regrOutput$p_il$importance,
                         p_nr = frst$regrOutput$p_nr$importance,
                         p_uk = frst$regrOutput$p_uk$importance)
imp_plt_dt$neat_names <- c("host response potential", "event order",
                           "res. transmission", "dens. transmission",
                           "background mort.", "infect. mort.",
                           "recovery rate",
                           "avg. reproduction",
                           "mutation prob.",
                           "carrying capacity", "carrying capacity SD",
                           "pre-disease time", 
                           "disease gens.",
                           "init. allele freq.", "init. infect.",
                           "post-disease time",
                           "adaptation pathway", "adaptive benefit",
                           "dominance", "trait SD", "adaptive cost")
imp_plt_long <- pivot_longer(imp_plt_dt, cols = 2:6, names_to = "outcome")
imp_plt_long$outcome <- recode(imp_plt_long$outcome,
                            p_er = "P(All ER)",
                            p_ex = "P(Ext)",
                            p_il = "P(IL)",
                            p_nr = "P(NR)",
                            p_uk = "Inconclusive")
imp_plt_long <- imp_plt_long %>% 
  mutate(fct_relevel(outcome, c("P(Ext)", "P(All ER)", "P(IL)",  "P(NR)", "Inconclusive")))

imp_plt_short <- imp_plt_long %>% # filtering down to top 8, that look appricably different from 0
  filter(neat_names %in% c("host response potential", 
                           "adaptation pathway", 
                           "res. transmission",
                           "dominance", 
                           "disease gens.",
                           "trait SD"))
tmp_plt_short <- imp_plt_short %>% filter(outcome %in% c("P(All ER)", "P(Ext)", "P(IL)"))
tmp_plt_short$outcome_f <- factor(tmp_plt_short$outcome, levels = c("P(Ext)", "P(All ER)", "P(IL)"))

fig1C <- ggplot(data = tmp_plt_short, aes(`value`, reorder(neat_names, `value`))) +
  geom_linerange(aes(xmin = 0, xmax = `value`)) +
  geom_point(aes(col = outcome), size = 2) + # alt: aes(size = log(IncNodePurity))
  labs(x = "Importance", y = NULL, col = "importance to event:") +
  scale_color_manual(values = c("#DB6341", "#ac1457", "#f1c4a2")) +
  facet_wrap(~outcome_f, nrow = 1) +
  labs(col = NULL) + 
  theme_bw(base_size = 9.5) + guides(col = "none") + 
  theme(legend.position = "bottom", 
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.key.spacing.y = unit(0, "pt"), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))  

fig1_h <- plot_grid(fig1A, fig1B, rel_widths = c(0.5, 1), labels = c("A", "B"))
fig1 <- plot_grid(fig1_h, fig1C, ncol = 1, rel_heights = c(1, 0.4), labels = c("", "C"))
# fig1

# save figure
pdf("figs/figure_plot/fig2_HP.pdf",height=160/25.4,width=170/25.4)
print(fig1)
dev.off()

# supp figures re: GSA
figS3 <- ggplot(data = imp_plt_long, aes(`value`, reorder(neat_names, `value`))) +
  geom_linerange(aes(xmin = 0, xmax = `value`)) +
  geom_point(aes(col = outcome), size = 3) + # alt: aes(size = log(IncNodePurity))
  labs(x = NULL, y = NULL, col = "Importance to event:") +
  scale_color_manual(values = c( "gray80", "#ac1457", "#DB6341",  "#f1c4a2", "gray40")) +
  facet_wrap(~outcome, nrow = 1) +
  theme_bw(base_size = 9.5) + guides(col = "none") + 
  theme(legend.position = "bottom", 
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.key.spacing.y = unit(0, "pt"), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))  

# save figure
pdf("figs/figure_plot/figS3_HP.pdf",height=140/25.4,width=170/25.4)
print(figS3)
dev.off()

#-----FIG2: COST/BENEFIT-----
sis_cb <- expand.grid("benefit" = c(100, 90, 80, 70, 0), # for percent change in benefit: 0, 50, 75, 100
                      "cost" = c("High cost\n(Mutant fecundity = 40%)", "Moderate cost\n(Mutant fecundity = 70%)", "No cost\n(Mutant fecundity = 100%)"), # 0.1 = some cost, 1 = no cost
                      "transmission" = c("Density only", "Density w/reservior", "Lasting disease shock"))
N <- dim(sis_cb)[1]
sis_cb$number <- 1:N
cost_ben <- merge(cb_dtf, sis_cb, by.x = "parm_number", by.y = "number", all = T) # %>% filter(cost > 1.8)

# get summaries of each probability to plot
ex_dat <- cost_ben %>% filter(extinct == 1) %>% 
  group_by(benefit, cost, transmission) %>% 
  summarize(`P(Ext)` = n()/2500)
# again, using "all ER" category
er_dat <- cost_ben %>% filter(extinct == 0 & abs(last_50-first_50-tot_50+1) < 0.8*tot_50 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 1) %>% 
  group_by(benefit, cost, transmission) %>% 
  summarize(`P(All ER)` = n()/2500, 
            `Time to\nrecovery (ER)` = mean((firstK95-120/3)/3), 
            upp_qnt = quantile((firstK95-120/3)/3, probs = 0.975),
            lwr_qnt = quantile((firstK95-120/3)/3, probs = 0.025))
# clean up time to recovery--recall that if there's no benefit, no evolution, so no ER to consider
er_dat$`Time to\nrecovery (ER)` <- ifelse(er_dat$benefit == 0, NA, er_dat$`Time to\nrecovery (ER)`)
er_dat$upp_qnt <- ifelse(er_dat$benefit == 0, NA, er_dat$upp_qnt)
er_dat$lwr_qnt <- ifelse(er_dat$benefit == 0, NA, er_dat$lwr_qnt)
dr_dat <- cost_ben %>% filter(extinct == 0 & abs(last_50-first_50-tot_50+1) < 0.8*tot_50 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 0 & final_inf_prev == 0) %>%    
  group_by(benefit, cost, transmission) %>% 
  summarize(`P(IL)` = n()/2500)

# merge data together
plt_dat <- merge(ex_dat, er_dat, all = T)
plt_dat <- merge(plt_dat, dr_dat, all = T) 
# clean up 
plt_dat$`P(Ext)` <- ifelse(is.na(plt_dat$`P(Ext)`), 0, plt_dat$`P(Ext)`)
plt_dat$`P(All ER)` <- ifelse(is.na(plt_dat$`P(All ER)`), 0, plt_dat$`P(All ER)`)
plt_dat$`P(IL)` <- ifelse(is.na(plt_dat$`P(IL)`), 0, plt_dat$`P(IL)`)
# rotate
cb_summ_long <- pivot_longer(plt_dat, cols= c(4:6, 9), names_to = "outcome", values_to = "value")
# add line range info
cb_summ_long$plus <- ifelse(cb_summ_long$outcome == "Time to\nrecovery (ER)", cb_summ_long$upp_qnt, cb_summ_long$value + intv)
cb_summ_long$less <- ifelse(cb_summ_long$outcome == "Time to\nrecovery (ER)", cb_summ_long$lwr_qnt, cb_summ_long$value - intv)
for (i in 1:length(cb_summ_long$value)) {
  if (cb_summ_long$plus[i] > 1 & cb_summ_long$outcome[i] != "Time to\nrecovery (ER)") cb_summ_long$plus[i] <- 1
  if (cb_summ_long$less[i] < 0 & cb_summ_long$outcome[i] != "Time to\nrecovery (ER)") cb_summ_long$less[i] <- 0
}
cb_summ_long$outcome_f <- factor(cb_summ_long$outcome, levels = c("P(Ext)", "P(All ER)", "P(IL)",
                                                            "Time to\nrecovery (ER)"))

fig2 <- ggplot(NULL, aes(x=as.factor(benefit), y=value, col = outcome_f, shape = as.factor(transmission))) + 
  geom_hline(data = cb_summ_long %>% filter(`benefit` == 0), aes(yintercept = value), col = "gray70", lty = "dashed") + 
  coord_cartesian(ylim = c(0, NA)) +
  geom_linerange(data = cb_summ_long, aes(ymin = less, ymax = plus), lwd = 1, position = position_dodge(width = 0.7)) + 
  geom_point(data = cb_summ_long, size = 2.25, fill = "white", position = position_dodge(width = 0.7)) + 
  scale_shape_manual(values = c(15, 16, 21)) + 
  scale_color_manual(values = c("#ac1457", "#DB6341", "#f1c4a2", "black")) +
  facet_grid(rows = vars(outcome_f), cols = vars(cost), scales = "free") + 
  guides(col = "none") + 
  theme_bw(base_size = 12) + 
  theme(legend.position = "bottom", 
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.key.spacing.y = unit(0, "pt")) + 
  labs(x = "Mutant allele benefit\n(percent change in rate)", shape = NULL, y = NULL)

# fig2

# # save figure
# png("figs/figure_plot/fig3_0604.png",height=160,width=170,res=400,units='mm')
# print(fig2)
# dev.off()

pdf("figs/figure_plot/fig3_HP.pdf",height=175/25.4,width=170/25.4)
print(fig2)
dev.off()

#-----FIG3: POP SIZES-----
# merge data
sis_pop_mut <- expand.grid("population size" = c(50, 100, 500),
                           "robustness" = "Mortality-blocking alleles", # M, G only bc they are strongest ER from past exploration
                           # "mutation rate" = c(0.0005, 0.005, 0.01),
                           "compartments" = "Recovery",
                           "transmission type" = c("Density only", "Density w/reservior", "Lasting disease shock"))
N <- dim(sis_pop_mut)[1]
sis_pop_mut$number <- 1:N
ex_risk <- merge(pm_dtf, sis_pop_mut, by.x = "parm_number", by.y = "number", all = T)

# get summary probabilities to plot
ex_dat <- ex_risk %>% filter(extinct == 1) %>% 
  group_by(parm_number) %>% 
  summarize(`P(Ext)` = n()/2500)
er_dat <- ex_risk %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) < 0.8*tot_50 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 1) %>% 
  group_by(parm_number) %>% 
  summarize(`P(All ER)` = n()/2500, 
            `Time to\nrecovery (ER)` = mean((firstK95-120/3)/3), 
            upp_qnt = quantile((firstK95-120/3)/3, probs = 0.975),
            lwr_qnt = quantile((firstK95-120/3)/3, probs = 0.025))
dr_dat <- ex_risk %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) < 0.8*tot_50 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 0 & final_inf_prev == 0) %>%    
  group_by(parm_number) %>% 
  summarize(`P(IL)` = n()/2500)

# merge
plt_dat <- merge(ex_dat, er_dat, all = T)
plt_dat <- merge(plt_dat, dr_dat, all = T)
tmp <- merge(plt_dat, sis_pop_mut, by.x = "parm_number", by.y = "number")
tmp$`P(IL)` <- ifelse(is.na(tmp$`P(IL)`), 0, tmp$`P(IL)`)

# rotate
pm_summ_long <- pivot_longer(tmp[, c(2:11)], cols= c(1:3, 6), names_to = "outcome", values_to = "value")
# get line range info
pm_summ_long$plus <- ifelse(pm_summ_long$outcome == "Time to\nrecovery (ER)", pm_summ_long$upp_qnt, pm_summ_long$value + intv)
pm_summ_long$less <- ifelse(pm_summ_long$outcome == "Time to\nrecovery (ER)", pm_summ_long$lwr_qnt, pm_summ_long$value - intv)
for (i in 1:length(pm_summ_long$value)) {
  if (pm_summ_long$plus[i] > 1 & pm_summ_long$outcome[i] != "Time to\nrecovery (ER)") pm_summ_long$plus[i] <- 1
  if (pm_summ_long$less[i] < 0 & pm_summ_long$outcome[i] != "Time to\nrecovery (ER)") pm_summ_long$less[i] <- 0
}
pm_summ_long$outcome_f <- factor(pm_summ_long$outcome, levels = c("P(Ext)", "P(All ER)", "P(IL)", 
                                                                  "Final M allele\nfrequency (ER)", "Final infection\nprevelence (ER)", "Time to\nrecovery (ER)"))

fig3 <- ggplot(NULL, aes(x=as.factor(`population size`), y=value, col = outcome_f, shape = as.factor(`transmission type`))) + 
  coord_cartesian(ylim = c(0, NA)) +
  geom_linerange(data = pm_summ_long, aes(ymin = less, ymax = plus), lwd = 1, position = position_dodge(width = 0.7)) + 
  geom_point(data = pm_summ_long, size = 2.25, fill = "white", position = position_dodge(width = 0.7)) + 
  scale_shape_manual(values = c(15, 16, 21)) + 
  scale_color_manual(values = c("#ac1457","#DB6341", "#f1c4a2", "black")) + 
  guides(col = "none") + 
  facet_grid(outcome_f ~ robustness, scales = "free_y") + 
  theme_bw(base_size = 12) + 
  theme(legend.position = "bottom", 
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.key.spacing.y = unit(0, "pt")) + 
  guides(shape = guide_legend(nrow = 2)) +
  labs(x = "Population size", shape = NULL, y = NULL)

# fig3

# # save figure
# png("figs/figure_plot/fig4_0604.png",height=160,width=85,res=400,units='mm')
# print(fig3)
# dev.off()
# # if 3x3 grid, width = 170 (this is pm_fig_A_0209 or pm_fig_B_0209, using mortality blocking (2))

pdf("figs/figure_plot/fig4_HP.pdf",height=175/25.4,width =85/25.4)
print(fig3)
dev.off()

#-----SUPP: DISEASE CYCLES-----
sis_dc <- expand.grid("disease cycles" = c(3:6), # aiming to span a bit of a range wihtout too many simulations
                      "robustness" = "Mortality-blocking", # M, B only
                      "event order" = c("Transmission,\nMortality\nRecovery", 
                                        "Mortality,\nTransmission\nRecovery", 
                                        "Recovery,\nMortality\nTransmission"),
                      "transmission" = c("Density only", "Density w/reservior", "Lasting disease shock"), 
                      "compartments" = "Recovery") # 1 = density, 2 = envrionemntal
N <- dim(sis_dc)[1]
sis_dc$parm_number <- 1:N
dc_dat <- merge(dc_dtf, sis_dc, by = "parm_number", all = T)

# get summary probabilities to plot
er_dat <- dc_dat %>% filter(extinct == 0 & abs(last_50-first_50-tot_50) < 0.8*tot_50 & tot_50 > 3 & at_K95 == 1 & r_allele_peak45 == 1) %>% 
  group_by(`parm_number`) %>% 
  summarize(`P(All ER)` = n()/2500, 
            `Time to\nrecovery (ER)` = mean((firstK95-(100/`disease cycles`))/`disease cycles`), 
            upp_qnt = quantile((firstK95-(100/`disease cycles`))/`disease cycles`, probs = 0.975),
            lwr_qnt = quantile((firstK95-(100/`disease cycles`))/`disease cycles`, probs = 0.025))
tmp <- merge(er_dat, sis_dc, by = "parm_number", all = T)
er_dat <- tmp; er_dat[is.na(er_dat)] <- 0
dr_dat <- dc_dat %>% filter(extinct == 0 & pop_drop50 == 1 & at_K95 == 1 & r_allele_peak45 == 0 & final_inf_prev == 0) %>%    
  group_by(parm_number) %>% 
  summarize(`P(IL)` = n()/2500)
tmp <- merge(dr_dat, sis_dc, by = "parm_number", all = T)
dr_dat <- tmp; dr_dat[is.na(dr_dat)] <- 0
ex_dat <- dc_dat %>% filter(extinct == 1) %>% 
  group_by(parm_number) %>% 
  summarize(`P(Ext)` = n()/2500)
tmp <- merge(ex_dat, sis_dc, by = "parm_number", all = T)
ex_dat <- tmp; ex_dat[is.na(ex_dat)] <- 0

# merge
tmp_dat <- merge(er_dat, ex_dat)
tmp_dat <- merge(tmp_dat, dr_dat)
# rotate
dc_all_long <- pivot_longer(tmp_dat[, 2:12], cols= c(6:7, 10:11), names_to = "outcome", values_to = "value")
# get line range info
dc_all_long$plus <- ifelse(dc_all_long$outcome == "Time to\nrecovery (ER)", dc_all_long$upp_qnt, dc_all_long$value + intv)
dc_all_long$less <- ifelse(dc_all_long$outcome == "Time to\nrecovery (ER)", dc_all_long$lwr_qnt, dc_all_long$value - intv)
for (i in 1:length(dc_all_long$value)) {
  if (dc_all_long$plus[i] > 1 & dc_all_long$outcome[i] != "Time to\nrecovery (ER)") dc_all_long$plus[i] <- 1
  if (dc_all_long$less[i] < 0 & dc_all_long$outcome[i] != "Time to\nrecovery (ER)") dc_all_long$less[i] <- 0
}
dc_all_long$outcome_f <- factor(dc_all_long$outcome, levels = c("P(Ext)", "P(All ER)", "P(IL)", 
                                                                "Final M allele\nfrequency (ER)", "Final infection\nprevelence (ER)", "Time to\nrecovery (ER)"))

figS4 <- ggplot(NULL, aes(x=as.factor(`disease cycles`), y=value, col = outcome_f, 
                                  shape=as.factor(transmission))) + 
  coord_cartesian(ylim = c(0, NA)) +
  geom_linerange(data = dc_all_long, aes(ymin = less, ymax = plus), lwd = 1, position = position_dodge(width = 0.7)) + 
  geom_point(data = dc_all_long, size = 2, fill = "white", position = position_dodge(width = 0.7)) + 
  scale_shape_manual(values = c(15, 16, 21)) + 
  facet_grid(row = vars(outcome_f), col = vars(`event order`), scales = "free_y") + 
  scale_color_manual(values = c("#ac1457","#DB6341", "#f1c4a2", "black")) + 
  guides(col = "none") + 
  theme_bw() + 
  theme_bw(base_size = 12) + 
  theme(legend.position = "bottom", 
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.key.spacing.y = unit(0, "pt")) + 
  guides(shape = guide_legend(nrow = 2)) +
  labs(x = "Disease cycles between reproduction", shape = NULL, y = NULL)

# figS4

# # save figure
# png("figs/figure_plot/SUPP_dc_0604.png",height=210,width=170,res=400,units='mm')
# print(figS4)
# dev.off()

pdf("figs/figure_plot/figS5_HP.pdf",height=175/25.4,width=170/25.4)
print(figS4)
dev.off()



