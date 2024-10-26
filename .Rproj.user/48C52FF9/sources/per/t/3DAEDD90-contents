library(tidyverse)
library(cowplot)
library(ggplot2)
library(viridis)

#-----FIG1: STRUCTURALS-----
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
# save figure
png("figures/figure_plot/structural_overview.png",height=155,width=170,res=400,units='mm')
print(fig1B)
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
plt_dat <- rbind(ex_dat, er_dat)
cb_summ_long <- pivot_longer(plt_dat, cols= 4:5, names_to = "outcome", values_to = "value")
cost_ben$firstK95 <- cost_ben$firstK95-100 # remove 100 year burn in period
cost_ben <- rename(cost_ben, 
                 `C - final\nM allele` = final_r_allele, 
                 `D - final\ninf prev` = final_inf_prev, 
                 `E - time\nto K` = firstK95)
cb_all_long <- pivot_longer(cost_ben %>% filter(extinct == 0 & pop_drop50 == 1 & at_K95 == 1 & r_allele_peak45 == 1) %>% 
                              select(c(8, 10, 18, 20:22)), 
                            cols = 1:3, names_to = "outcome", values_to = "value")
fig2 <- ggplot(NULL, aes(x=as.factor(benefit), y=value, col = as.factor(cost))) + 
  geom_boxplot(data = cb_all_long) + 
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
fig3 <- ggplot(NULL, aes(x=as.factor(`population size`), y=value, col=as.factor(`mutation rate`))) + 
  geom_boxplot(data = pm_all_long) + 
  # coord_cartesian(ylim = c(0, NA)) + 
  geom_point(data = pm_summ_long, size = 3) + 
  scale_color_manual(values = c("#ac1457","#DB6341", "#f1c4a2")) + 
  facet_grid(row = vars(outcome), col = vars(`robustness`), scales = "free_y") + 
  theme_bw() + theme(text = element_text(size = 12), legend.position = "bottom") + 
  labs(x = "Population size", col = "Mutation rate", y = NULL)
# save figure
png("figures/figure_plot/ext_risk.png",height=115,width=170,res=400,units='mm')
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
