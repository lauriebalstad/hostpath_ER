library(tidyverse)
library(cowplot)
library(ggplot2)
library(viridis)

#-----FIG1: STRUCTURALS-----
# load in data
str_dat <- readRDS("figures/figure_data/all_ts_dat.rds")
# rotate data
all_ts_long <- pivot_longer(all_ts_dat %>% select(c(1, 7:10)), cols = c(2, 4, 5), names_to = "outcome", values_to = "value")
colors <- c("Extinction" = alpha("black", 0.8), "Persist. ER" = alpha("#ac1457",0.8),
            "Temp. ER" = "#f1c4a2", "DR" = alpha("#DB6341", 0.8))
fig1A <- ggplot(all_ts_long, aes(ts, value)) + # 7, 8, 9, 10 gives exinction
  geom_line(aes(group = trial), col = alpha("black", 0.1), lwd = 0.75) +
  geom_line(data = all_ts_long %>% filter(trial == 3), aes(color = "Extinction"), lwd = 0.85) + # extinction
  geom_line(data = all_ts_long %>% filter(trial == 11), aes(color = "Persist. ER"), lwd = 0.85) + # ER
  geom_line(data = all_ts_long %>% filter(trial == 12), aes(color = "Temp. ER"), lwd = 0.85) + # back selection
  geom_line(data = all_ts_long %>% filter(trial == 10), aes(color = "DR"), lwd = 0.85) + # DR
  facet_wrap(~outcome, scale = "free_y", ncol = 1) + theme_bw() + 
  scale_color_manual(values = colors) + 
  theme(text = element_text(size = 12), legend.position = "bottom") + 
  guides(color = guide_legend(nrow = 2)) +
  labs(x = "Standardized generation", y = NULL, col = NULL)

# load in data
str_dat <- readRDS("figures/figure_data/str_dat.rds")
# get summaries to plot
prob_dat <- str_dat %>% group_by(compartments, `transmission type`, robustness) %>% 
  summarize(`B - P(extinct)` = sum(extinct)/n(), #,
            `A - P(decline)` = sum(pop_drop50)/n()) # note for frequency dependent, force of infection is much much lower... think about how to make these more equalivalent
er_dat <- str_dat %>% filter(extinct == 0 & pop_drop50 == 1 & at_K95 == 1 & r_allele_peak45 == 1) %>% 
  group_by(compartments, `transmission type`, robustness) %>% 
  summarize(`D - P(rescue:\nevolutionary)` = n()/500)
dr_dat <- str_dat %>% filter(extinct == 0 & pop_drop50 == 1 & at_K95 == 1 & r_allele_peak45 == 0) %>% 
  group_by(compartments, `transmission type`, robustness) %>% 
  summarize(`C - P(rescue:\ndemographic)` = n()/500)
ma_dat <- str_dat %>% filter(extinct == 0 & pop_drop50 == 0 & at_K95 == 1 & r_allele_peak45 == 1) %>% 
  group_by(compartments, `transmission type`, robustness) %>% 
  summarize(`E - P(micro-\nadaptation)` = n()/500)
# merge
tmp <- merge(prob_dat, er_dat, by = c("compartments", "transmission type", "robustness"), all = T)
tmp <- merge(tmp, ma_dat, by = c("compartments", "transmission type", "robustness"), all = T)
plot_dat <- merge(tmp, dr_dat, by = c("compartments", "transmission type", "robustness"), all = T)
# rename
plot_dat$compartments <- recode(plot_dat$compartments, "1" = "SIX", "2" = "SIS", "3" = "SIR")
plot_dat$`transmission type` <- recode(plot_dat$`transmission type`, "1" = "density", "2" = "environmental")
plot_dat$robustness <- recode(plot_dat$robustness, "1" = "\u03bc", "2" = "\u03b2", "3" = "\u03d2", "4" = "N")
# convert to long
plot_long <- pivot_longer(plot_dat, cols = 4:8, names_to = "outcome", values_to = "value")
plot_long$value <- ifelse(is.na(plot_long$value), 0, plot_long$value)
# points plot
fig1B <- ggplot(plot_long, aes(robustness, value, col = `transmission type`)) + 
  geom_hline(data = plot_long %>% filter(robustness == "N"), aes(yintercept = value), col = "gray70", lty = "dashed") + 
  geom_point(size = 3) + scale_color_manual(values = c("#ac1457", "#f1c4a2")) +
  facet_grid(rows = vars(outcome), cols = vars(compartments)) + 
  labs(x = "Robustness type", y = "Probability", col = "Transmission type") + 
  theme_bw()  + 
  theme(text = element_text(size = 12), legend.position = "bottom") 

# combine figure 1
fig1 <- plot_grid(fig1A, fig1B, rel_widths= c(0.5, 1), nrow = 1)

# save figure
png("figures/figure_plot/structural_overview.png",height=155,width=170,res=400,units='mm')
print(fig1)
dev.off()

#-----FIG2: COST/BENEFIT-----
# load in data
cost_ben <- readRDS("figures/figure_data/cost_ben_dat.rds")
# get summaries to plot -- P(ER) + T_K + max/final R/inf (6) for the ER cases only
ex_dat <- cost_ben %>% filter(extinct == 1) %>% 
  group_by(benefit, cost, case) %>% 
  summarize(`A - P(extinct)` = n()/500)
er_dat <- cost_ben %>% filter(extinct == 0 & pop_drop50 == 1 & at_K95 == 1 & r_allele_peak45 == 1) %>% 
  group_by(benefit, cost, case) %>% 
  summarize(`B - P(rescue:\nevolutionary)` = n()/500)
# dr_dat <- cost_ben %>% filter(extinct == 0 & pop_drop50 == 1 & at_K95 == 1 & r_allele_peak45 == 0) %>% 
#   group_by(benefit, cost, case) %>% 
#   summarize(`C - P(rescue:\ndemographic)` = n()/500)
plt_dat <- rbind(ex_dat, er_dat, dr_dat)
cb_summ_long <- pivot_longer(plt_dat, cols= 4:6, names_to = "outcome", values_to = "value")
cost_ben$firstK95 <- cost_ben$firstK95-100 # remove 100 year burn in period
cost_ben <- rename(cost_ben, 
                 `C - final\nM allele` = final_r_allele, 
                 `D - final\ninf prev` = final_inf_prev, 
                 `E - time\nto K` = firstK95)
cb_all_long <- pivot_longer(cost_ben %>% filter(extinct == 0 & pop_drop50 == 1 & at_K95 == 1 & r_allele_peak45 == 1) %>% 
                              select(c(8, 10, 18, 20:22)), 
                            cols = 1:3, names_to = "outcome", values_to = "value")
# remove the no-evo cases from box plots
cb_mu_beta <- cb_all_long %>% filter(case %in% c("\u03bc", "\u03b2")) %>% filter(benefit < 1)
cb_gamma <- cb_all_long %>% filter(case == "\u03d2") %>% filter(benefit > 0)
cb_no_evo <- rbind(cb_mu_beta, cb_gamma)
fig2 <- ggplot(NULL, aes(x=as.factor(benefit), y=value, col = as.factor(cost))) + 
  geom_boxplot(data = cb_no_evo) + 
  geom_point(data = cb_summ_long, size = 3) + 
  scale_color_manual(values = c("#ac1457", "#f1c4a2")) + 
  facet_grid(rows = vars(outcome), cols = vars(case), scales = "free_y") + 
  theme_bw() + theme(text = element_text(size = 12), legend.position = "bottom") + 
  labs(x = "Allele strength", col = "Allele fecundity", y = NULL)
# save figure
png("figures/figure_plot/cost_benefit.png",height=155,width=170,res=400,units='mm')
print(fig2)
dev.off()

#-----FIG3: POP SIZES-----
# load in data
ex_risk <- readRDS("figures/figure_data/ext_risk_dat.rds")
# rename
ex_risk$robustness <- recode(ex_risk$robustness, "1" = "\u03bc", "2" = "\u03b2", "3" = "\u03d2", "4" = "N")
# get summaries to plot -- P(ER) + T_K + max/final R/inf (6) for the ER cases only
ex_dat <- ex_risk %>% filter(extinct == 1) %>% 
  group_by(`population size`, `robustness`, `mutation rate`) %>% 
  summarize(`A - P(extinct)` = n()/500)
er_dat <- ex_risk %>% filter(extinct == 0 & pop_drop50 == 1 & at_K95 == 1 & r_allele_peak45 == 1) %>% 
  group_by(`population size`, `robustness`, `mutation rate`) %>% 
  summarize(`B - P(rescue:\nevolutionary)` = n()/500)
plt_dat <- rbind(ex_dat, er_dat)
pm_summ_long <- pivot_longer(plt_dat, cols= 4:5, names_to = "outcome", values_to = "value")
ex_risk$firstK95 <- ex_risk$firstK95-100 # remove 100 year burn in period
ex_risk <- rename(ex_risk, 
                   `C - time\nto K` = firstK95)
pm_all_long <- pivot_longer(ex_risk %>% filter(extinct == 0 & pop_drop50 == 1 & at_K95 == 1 & r_allele_peak45 == 1) %>% 
                              select(c(18, 20:22)), 
                            cols = 1, names_to = "outcome", values_to = "value")
figS2 <- ggplot(NULL, aes(x=as.factor(`population size`), y=value, col=as.factor(`mutation rate`))) + 
  geom_boxplot(data = pm_all_long) + 
  # coord_cartesian(ylim = c(0, NA)) + 
  geom_point(data = pm_summ_long, size = 3) + 
  scale_color_manual(values = c("#ac1457","#DB6341", "#f1c4a2")) + 
  facet_grid(row = vars(outcome), col = vars(`robustness`), scales = "free_y") + 
  theme_bw() + theme(text = element_text(size = 12), legend.position = "bottom") + 
  labs(x = "Population size", col = "Mutation rate", y = NULL)
# save figure
png("figures/figure_plot/supp_ext_risk.png",height=115,width=170,res=400,units='mm')
print(figS2)
dev.off()

fig3 <- ggplot(NULL, aes(x=as.factor(`population size`), y=value, col=as.factor(`mutation rate`))) + 
  geom_boxplot(data = pm_all_long %>% filter(robustness == "\u03d2")) + 
  # coord_cartesian(ylim = c(0, NA)) + 
  geom_point(data = pm_summ_long %>% filter(robustness == "\u03d2"), size = 3) + 
  scale_color_manual(values = c("#ac1457","#DB6341", "#f1c4a2")) + 
  facet_grid(row = vars(outcome), col = vars(`robustness`), scales = "free_y") + 
  theme_bw() + theme(text = element_text(size = 12), legend.position = "bottom") + 
  labs(x = "Population size", col = "Mutation rate", y = NULL)
# save figure
png("figures/figure_plot/ext_risk.png",height=130,width=85,res=400,units='mm')
print(fig3)
dev.off()

#-----SUPP: DISEASE CYCLES-----
# load in data
dc_dat <- readRDS("figures/figure_data/cyc_ord_dat.rds")
# rename
dc_dat$`event order` <- recode(dc_dat$`event order`, "2" = "\u03b2, \u03d2, \u03bc", "4" = "\u03bc, \u03b2, \u03d2", "6" = "\u03d2, \u03bc, \u03b2")
dc_dat$robustness <- recode(dc_dat$robustness, "1" = "\u03bc", "2" = "\u03b2", "3" = "\u03d2", "4" = "N")
# get summaries to plot -- P(ER) + T_K + max/final R/inf (6) for the ER cases only --> note more generous def of pop drop
er_dat <- dc_dat %>% filter(extinct == 0 & pop_drop80 == 1 & at_K95 == 1 & r_allele_peak15 == 1) %>% 
  group_by(`disease cycles`, `robustness`, `event order`, `transmission`) %>% 
  summarize(`B - P(rescue:\nevolutionary)` = n()/500, 
            `final R allele` = mean(final_r_allele), 
            `max R allele` = mean(max_r_allele),
            `final inf prev` = mean(final_inf_prev),
            `max inf prev` = mean(max_inf_prev), 
            `time to K` = mean(firstK95))
ex_dat <- dc_dat %>% filter(extinct == 1) %>% group_by(`disease cycles`, `robustness`, `event order`, `transmission`) %>% 
  summarize(`A - P(extinct)` = n()/500)
dc_er_long <- pivot_longer(er_dat, cols= c(5), names_to = "outcome", values_to = "value")
dc_ex_long <- pivot_longer(ex_dat, cols= c(5), names_to = "outcome", values_to = "value")
dc_dat <- rename(dc_dat, 
                   `C - final\nM allele` = final_r_allele, 
                   `D - final\ninf prev` = final_inf_prev, 
                   `E - time\nto K` = firstK95)
dc_all_long <- pivot_longer(dc_dat %>% filter(extinct == 0 & pop_drop80 == 1 & at_K95 == 1 & r_allele_peak15 == 1) %>% 
                              select(c(8, 10, 18, 20:23)), 
                            cols = 1:3, names_to = "outcome", values_to = "value")
fig_S1 <- ggplot(NULL, aes(x=as.factor(`disease cycles`), y=value, col=as.factor(robustness))) + 
  geom_boxplot(data = dc_all_long) + coord_cartesian(ylim = c(0, NA)) + 
  geom_point(data = dc_er_long, size = 3) + 
  geom_point(data = dc_ex_long, size = 3) + 
  scale_color_manual(values = c("#ac1457", "#f1c4a2", "gray70")) + 
  facet_grid(row = vars(outcome), col = vars(`event order`), scales = "free_y") + 
  theme_bw() + theme(text = element_text(size = 12), legend.position = "bottom") + 
  labs(x = "Disease cycles between reproduction", col = "Robustness type", y = NULL)
# save figure
png("figures/figure_plot/supp_disease_eco.png",height=155,width=170,res=400,units='mm')
print(fig_S1)
dev.off()


#-----NOT USED: INF FORCE-----
# load in data
inf_force <- readRDS("figures/figure_data/dis_force_dat.rds")
# rename
inf_force$`transmission type` <- recode(inf_force$`transmission type`, "1" = "density", "2" = "environmental")
# get summaries to plot -- P(ER) + T_K + max/final R/inf (6) for the ER cases only
ex_dat <- inf_force %>% filter(extinct == 1) %>% 
  group_by(`transmission type`, `base transmission rate`, `mortality rate`) %>% 
  summarize(`A - P(extinct)` = n()/500)
er_dat <- inf_force %>% filter(extinct == 0 & pop_drop20 == 1 & at_K95 == 1 & r_allele_peak75 == 1) %>% 
  group_by(`transmission type`, `base transmission rate`, `mortality rate`) %>% 
  summarize(`B - P(rescue:\nevolutionary)` = n()/500)
dr_dat <- inf_force %>% filter(extinct == 0 & pop_drop20 == 1 & at_K95 == 1 & r_allele_peak75 == 0) %>% 
  group_by(`transmission type`, `base transmission rate`, `mortality rate`) %>% 
  summarize(`C - P(rescue:\ndemographic)` = n()/500)
plt_dat <- rbind(ex_dat, er_dat, dr_dat)
plt_dat_75 <- data.frame(rep("density", 2), c(1, 2), c(0.75, 0.75), rep(0, 2), rep(0, 2), rep(0, 2), rep(NA, 2))
colnames(plt_dat_75) <- colnames(plt_dat)
plt_dat <- rbind(plt_dat, plt_dat_75)
pm_summ_long <- pivot_longer(plt_dat, cols= 4:6, names_to = "outcome", values_to = "value")
pm_summ_long$robustness <- rep("\u03d2")
fig4 <- ggplot(NULL, aes(x=as.factor(`mortality rate`), y=value, col=as.factor(`base transmission rate`))) + 
  geom_point(data = pm_summ_long, size = 3) + 
  coord_cartesian(ylim = c(0, 1)) + 
  scale_color_manual(values = c("#ac1457", "#f1c4a2")) + 
  facet_grid(row = vars(outcome), col = vars(`robustness`), scales = "free_y") + 
  theme_bw() + theme(text = element_text(size = 12), legend.position = "bottom") + 
  labs(x = "Mortality rate", col = "Transmission type", y = NULL)
# save figure
png("figures/figure_plot/infection_force.png",height=115,width=170,res=400,units='mm')
print(fig4)
dev.off()
