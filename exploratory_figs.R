library(tidyverse)
library(cowplot)
library(ggplot2)
library(viridis)

base_vect <- c(1, 2, 1, # SIX, BGM, density
               15, 10, # timing things
               1, 1, 1, 0.05, # 6-9: transmission
               0.05, 0.05, 0.05, # background mort + 16
               1, 1, 1, # 13-15 + 17: disease mort
               0.01, 0.02, # mort sd
               0.05, 0.05, 0.05, 0.01, # 18-21: recovery, unused bc SIX (1)
               1.5, 1.75, 2, 0.2, 0.01, # reproduction & mutation
               100, 4, # carrying capacity things 
               1) # nubmer of disease cycles per gen -- tyring 2
# nb: in base case, no benefit of allele, but cost is fixed

# so.... fig 1 will be comparision of structure across base model
cases <- expand.grid("transmission type" = 1:2, 
                     "compartments" = 1:3, # note SIR doesn't seem to have a big effect
                     "robustness" = 1:4) # 1 = mortality, 2 = transmission, 3 = recovery, 4 = demographic rescue only
# remove 1/3 compartment/robustness combo
cases <- cases[c(1:12, 15:24), ]
sim_dat <- NULL # this is where everything will get stored

# doing parallel
cl <- parallel::makeCluster(detectCores())
doParallel::registerDoParallel(cl)

for (case in 1:dim(cases)[1]){
  
  print(case) # for error id and to check in
  
  case_vect <- base_vect
  # if transmission type is frequency (2), need to increase baseline transmission a bit
  if (cases$`transmission type`[case] == 2) {
    case_vect[6:8] <- 2 # transmission increased
  }
  if (cases$robustness[case] == 1) {
    case_vect[13] <- 0
    case_vect[14] <- 0.5
  } # 1 = mortality
  if (cases$robustness[case] == 2) {
    case_vect[6] <- 0*case_vect[8] # since will be scaled by WW rate
    case_vect[7] <- 0.5*case_vect[8]
  } # 1 = transmission
  if (cases$robustness[case] == 3) {
    case_vect[13] <- 1
    case_vect[14] <- 0.5
  } # 1 = recovery, nb: recovery higher
  
  if (cases$robustness[case] == 4) {
    
    case_vect <- base_vect # no robustness
    
    # demographics only...
    case_vect[22:24] <- 2 # no cost 
    
    # and make sure that if frequency dependent, things are set up correctly
    if (cases$`transmission type`[case] == 2) {
      case_vect[6:8] <- 3 # WW transmission increased
    }
    
  } # make sure that everything is reset...
  
  # modify base parameters
  case_vect[3] <- cases$`transmission type`[case] # since numbers have meaning
  case_vect[1] <- cases$compartments[case] # diddo
  
  # run simulation
  case_dat <- run_gens(case_vect, 70, 500)
  
  # save connect to case info
  case_dat$`transmission type` <- rep(cases$`transmission type`[case])
  case_dat$compartments <- rep(cases$compartments[case])
  case_dat$robustness <- rep(cases$robustness[case])
  
  # then append to sim_dat
  sim_dat <- rbind(sim_dat, case_dat)
    
}

# end parallel cluster
stopCluster(cl)

fig1_dat <- sim_dat # G = 0.05
fig1_dat_highG <- sim_dat # G = 0.1

# save things!
saveRDS(fig1_dat, file = "figure_data/fig1_dat.rds")
saveRDS(fig1_dat_highG, file = "figure_data/fig1_dat_highG.rds")

# get summaries to plot
prob_dat <- sim_dat %>% group_by(compartments, `transmission type`, robustness) %>% 
  summarize(`P(extinct)` = sum(extinct)/n(), 
            `P(decline)` = sum(pop_drop)/n(), # note for frequency dependent, force of infection is much much lower... think about how to make these more equalivalent
            `P(population at K)` = sum(at_K)/n(), # did the population recover?
            `P(R allele > 0.4)` = sum(r_allele_peak)/n(), 
            `P(lost disease)` = length(which(final_inf_prev == 0))/n())
er_dat <- sim_dat %>% filter(extinct == 0 & pop_drop == 1 & at_K == 1 & r_allele_peak == 1) %>% 
  group_by(compartments, `transmission type`, robustness) %>% 
  summarize(`P(rescue: evolutionary)` = n()/500)
dr_dat <- sim_dat %>% filter(extinct == 0 & pop_drop == 1 & at_K == 1 & r_allele_peak == 0) %>% 
  group_by(compartments, `transmission type`, robustness) %>% 
  summarize(`P(rescue: demographic)` = n()/500)
# merge
tmp <- merge(prob_dat, er_dat, by = c("compartments", "transmission type", "robustness"), all = T)
plot_dat <- merge(tmp, dr_dat, by = c("compartments", "transmission type", "robustness"), all = T)
# rename
plot_dat$compartments <- recode(plot_dat$compartments, "1" = "SIX", "2" = "SIS", "3" = "SIR")
plot_dat$`transmission type` <- recode(plot_dat$`transmission type`, "1" = "density", "2" = "frequency")
plot_dat$robustness <- recode(plot_dat$robustness, "1" = "M", "2" = "B", "3" = "G", "4" = "N")
# convert to long
plot_long <- pivot_longer(plot_dat, cols = 4:10, names_to = "outcome", values_to = "value")
plot_long$value <- ifelse(is.na(plot_long$value), 0, plot_long$value)
# points plot -- also show varience? or issue with number of simulations?
ggplot(plot_long, aes(robustness, value, col = `transmission type`)) + 
  geom_hline(data = plot_long %>% filter(robustness == "N"), aes(yintercept = value), col = "gray70", lty = "dashed") + 
  geom_point(size = 2) + scale_color_manual(values = c("#ac1457", "#f1c4a2")) +
  facet_grid(rows = vars(compartments), cols = vars(outcome)) + 
  labs(x = "evolutionary benefit", y = "probability") + 
  theme_bw() # + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# for box plots, need full data set!
tmp <- sim_dat %>% filter(extinct == 0 & pop_drop == 1 & at_K == 1 & r_allele_peak == 1 & robustness != 4) 
sim_long <- pivot_longer(tmp, cols = 1:14, names_to = "outcome", values_to = "value") %>% 
  filter(outcome %in% c("final_r_allele", "max_r_allele", 
                        "final_inf_prev", "max_inf_prev", 
                        "final_pop_size", "first_K")) # since only ER cases, don't really care about pop outcomes? should get time first at K
sim_long$compartments <- recode(sim_long$compartments, "1" = "SIX", "2" = "SIS", "3" = "SIR")
sim_long$`transmission type` <- recode(sim_long$`transmission type`, "1" = "density", "2" = "frequency")
sim_long$robustness <- recode(sim_long$robustness, "1" = "M", "2" = "B", "3" = "G", "4" = "N")
sim_long$outcome <- recode(sim_long$outcome, "final_r_allele" = "T_95 freq: R allele", "max_r_allele" = "Max. freq: R allele", 
                                             "final_inf_prev" = "T_95 freq: Infection", "max_inf_prev" = "Max. freq: Infection", 
                                             "final_pop_size" = "T_95 population size", "first_K" = "T_K")
six_plt <- ggplot(sim_long %>% filter(compartments == "SIX"), aes(x = robustness, y = value, col = `transmission type`)) + 
  # geom_hline(data = plot_long %>% filter(robustness == "N"), aes(yintercept = median(value)), col = "gray70", lty = "dashed") + 
  geom_boxplot() + scale_color_manual(values = c("#ac1457", "#f1c4a2")) +
  facet_wrap(~outcome, scales = "free_y", nrow = 3) + 
  labs(x = "evolutionary benefit", y = NULL) + 
  theme_bw() + theme(text = element_text(size = 12))
sis_plt <- ggplot(sim_long %>% filter(compartments == "SIS"), aes(x = robustness, y = value, col = `transmission type`)) + 
  # geom_hline(data = plot_long %>% filter(robustness == "N"), aes(yintercept = median(value)), col = "gray70", lty = "dashed") + 
  geom_boxplot() + scale_color_manual(values = c("#ac1457", "#f1c4a2")) +
  facet_wrap(~outcome, scales = "free_y", nrow = 6) + 
  labs(x = "evolutionary benefit", y = NULL) + 
  theme_bw() + theme(text = element_text(size = 12), legend.position = "bottom")
sir_plt <- ggplot(sim_long %>% filter(compartments == "SIR"), aes(x = robustness, y = value, col = `transmission type`)) + 
  # geom_hline(data = plot_long %>% filter(robustness == "N"), aes(yintercept = median(value)), col = "gray70", lty = "dashed") + 
  geom_boxplot() + scale_color_manual(values = c("#ac1457", "#f1c4a2")) +
  facet_wrap(~outcome, scales = "free_y", nrow = 3) + 
  labs(x = "evolutionary benefit", y = NULL) + 
  theme_bw() + theme(text = element_text(size = 12))
leg_plt <- get_legend(six_plt + theme(legend.box.margin = margin(0, 0, 0, 12)))
plot_grid(six_plt + theme(legend.position="none"), 
          sis_plt + theme(legend.position="none"), 
          leg_plt, nrow = 1, rel_widths = c(1, 1, 0.35))

# comparing different cost/benefits
# so.... fig 1 will be comparision of structure across base model
sis_cb <- expand.grid("benefit" = c(0, 0.5, 1, 1.5), # 0 = total benefit, 1.5 = some cost for B & M, but opposite for G
                        "cost" = c(1.25, 2), # 1.25 = some cost, 2 = no cost
                        "case" = c("B", "M", "G")) 
cb_dat <- NULL # this is where everything will get stored

# doing parallel
cl <- parallel::makeCluster(detectCores())
doParallel::registerDoParallel(cl)

for (cb in 1:dim(sis_cb)[1]){
  
  print(cb) # for error id and to check in
  
  cb_vect <- base_vect
  cb_vect[1] <- 2 # SIS, with density
  # benefit update
  if (sis_cb$case[cb] == "B") {
    cb_vect[6] <- sis_cb$benefit[cb] # sis_cb$benefit[cb]*cb_vect[8]
    cb_vect[7] <- (cb_vect[6]+cb_vect[8])/2}
  if (sis_cb$case[cb] == "M") { 
    cb_vect[13] <- sis_cb$benefit[cb] # sis_cb$benefit[cb]*cb_vect[15]
    cb_vect[14] <- (cb_vect[13]+cb_vect[15])/2}
  if (sis_cb$case[cb] == "G") { 
    cb_vect[18] <- sis_cb$benefit[cb] # sis_cb$benefit[cb]*cb_vect[20]
    cb_vect[19] <- (cb_vect[18]+cb_vect[20])/2}
  # note for G, the meaning is opposite
  
  # cost update -- not dependent on the case
  cb_vect[22] <- sis_cb$cost[cb]
  cb_vect[23] <- (cb_vect[22]+cb_vect[24])/2
  
  # run simulation
  case_dat <- run_gens(cb_vect, 70, 500)
  
  # save connect to case info
  case_dat$benefit <- rep(sis_cb$benefit[cb])
  case_dat$cost <- rep(sis_cb$cost[cb])
  case_dat$case <- rep(sis_cb$case[cb])

  # then append to sim_dat
  cb_dat <- rbind(cb_dat, case_dat)
  
}

# end parallel cluster
stopCluster(cl)

fig2_dat_highG <- cb_dat # G = 0.1
fig2_dat <- cb_dat # G = 0.05

# save things!
saveRDS(fig2_dat, file = "figure_data/fig2_dat.rds")
saveRDS(fig2_dat_highG, file = "figure_data/fig2_dat_highG.rds")

# get summaries to plot -- P(ER) + T_K + max/final R/inf (6) for the ER cases only
er_dat <- cb_dat %>% filter(extinct == 0 & pop_drop == 1 & at_K == 1 & r_allele_peak == 1) %>% 
  group_by(benefit, cost, case) %>% 
  summarize(`P(rescue: evolutionary)` = n()/500, 
            `final R allele` = mean(final_r_allele), 
            `max R allele` = mean(max_r_allele),
            `final inf prev` = mean(final_inf_prev),
            `max inf prev` = mean(max_inf_prev), 
            `time to K` = mean(first_K))
cb_summ_long <- pivot_longer(er_dat, cols= 4:9, names_to = "outcome", values_to = "value")
cb_all_long <- pivot_longer(cb_dat %>% filter(extinct == 0 & pop_drop == 1 & at_K == 1 & r_allele_peak == 1) %>% select(c(5, 11, 7:9, 15:17)), 
                            cols = 1:5, names_to = "outcome", values_to = "value")
ggplot(NULL, aes(x=as.factor(benefit), y=value, col = as.factor(cost))) + 
  geom_boxplot(data = cb_all_long) + 
  geom_point(data = cb_summ_long %>% filter(outcome == "P(rescue: evolutionary)"), size = 3) + scale_color_manual(values = c("#ac1457", "#f1c4a2")) + 
  facet_grid(rows = vars(outcome), cols = vars(case), scales = "free_y") + 
  theme_bw() + theme(text = element_text(size = 12)) + labs(x = "R allele strength", col = "R allele\nfecundity", y = NULL)

# then will do time series plot (note need a slightly different funciton to create these)

# other variable to check in the number of disease cycles between reproductions
sis_dc <- expand.grid("disease cycles" = c(1, 2, 3), # 0 = total benefit, 1.5 = some cost for B & M, but opposite for G
                      "case" = c("B", "M", "G")) 
dc_dat <- NULL # this is where everything will get stored

# doing parallel
cl <- parallel::makeCluster(detectCores())
doParallel::registerDoParallel(cl)

for (dc in 1:dim(sis_dc)[1]){
  
  print(dc) # for error id and to check in
  
  dc_vect <- base_vect
  dc_vect[1] <- 2 # SIS, with density
  # need to redefine each component based on number of cycles... talk to marissa about how to do this
  dc_vect[6:8] <- base_vect[6:8]/sis_dc$`disease cycles`[dc]
  dc_vect[13:15] <- base_vect[13:15]/sis_dc$`disease cycles`[dc]
  dc_vect[18:20] <- base_vect[18:20]/sis_dc$`disease cycles`[dc]
  # benefit update
  if (sis_dc$case[dc] == "B") {
    dc_vect[6] <- 0 
    dc_vect[7] <- (dc_vect[6]+dc_vect[8])/2}
  if (sis_dc$case[dc] == "M") { 
    dc_vect[13] <- 0 
    dc_vect[14] <- (dc_vect[13]+dc_vect[15])/2}
  if (sis_dc$case[dc] == "G") { 
    dc_vect[18] <- 1/sis_dc$`disease cycles`[dc]
    dc_vect[19] <- (dc_vect[18]+dc_vect[20])/2}
  # note for G, the meaning is opposite
  
  # cost update -- not dependent on the case
  dc_vect[22] <- 2
  dc_vect[23] <- (dc_vect[22]+dc_vect[24])/2
  
  # disease cycle update
  dc_vect[29] <- sis_dc$`disease cycles`[dc]
  
  # run simulation
  case_dat <- run_gens(dc_vect, 70, 500)
  
  # save connect to case info
  case_dat$`disease cycles` <- rep(sis_dc$`disease cycles`[dc])
  case_dat$case <- rep(sis_dc$case[dc])
  
  # then append to sim_dat
  dc_dat <- rbind(dc_dat, case_dat)
  
}

# end parallel cluster
stopCluster(cl)

fig3_dat_highG <- dc_dat # G = 0.1
fig3_dat <- dc_dat # G = 0.05

# save things!
saveRDS(fig3_dat, file = "figure_data/fig3_dat.rds")
saveRDS(fig3_dat_highG, file = "figure_data/fig3_dat_highG.rds")

# get summaries to plot -- P(ER) + T_K + max/final R/inf (6) for the ER cases only
er_dat <- dc_dat %>% filter(extinct == 0 & pop_drop == 1 & at_K == 1 & r_allele_peak == 1) %>% 
  group_by(`disease cycles`, case) %>% 
  summarize(`P(rescue: evolutionary)` = n()/500, 
            `final R allele` = mean(final_r_allele), 
            `max R allele` = mean(max_r_allele),
            `final inf prev` = mean(final_inf_prev),
            `max inf prev` = mean(max_inf_prev), 
            `time to K` = mean(first_K))
dc_summ_long <- pivot_longer(er_dat, cols= 3:8, names_to = "outcome", values_to = "value")
dc_all_long <- pivot_longer(dc_dat %>% filter(extinct == 0 & pop_drop == 1 & at_K == 1 & r_allele_peak == 1) %>% select(c(5, 11, 7:9, 15:16)), 
                            cols = 1:5, names_to = "outcome", values_to = "value")
ggplot(NULL, aes(x=as.factor(`disease cycles`), y=value, col=case)) + 
  geom_boxplot(data = dc_all_long) + 
  geom_point(data = dc_summ_long %>% filter(outcome == "P(rescue: evolutionary)"), size = 3) + scale_color_manual(values = c("#ac1457", "#f1c4a2", "gray70")) + 
  facet_grid(row = vars(outcome), col = vars(case), scales = "free_y") + 
  theme_bw() + theme(text = element_text(size = 12)) + labs(x = "Disease cycles between reproduction", col = NULL, y = NULL)

