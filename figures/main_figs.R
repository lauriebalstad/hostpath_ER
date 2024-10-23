library(tidyverse)
library(cowplot)
library(ggplot2)
library(viridis)

base_vect <- c(2, 2, 1, # SIS, BGM, density
               100, 0, # timing things -- no wait between mutation and non-mutation
               1, 1, 1, 0.05, # 6-9: transmission
               0.05, 0.05, 0.05, # background mort + 16
               1, 1, 1, # 13-15 + 17: disease mort
               0.01, 0.02, # mort sd
               0.05, 0.05, 0.05, 0.01, # 18-21: recovery, unused bc SIX (1)
               1.5, 1.75, 2, 0.2, 0.005, # reproduction & mutation -- try smaller number
               100, 4, # carrying capacity things 
               1) # nubmer of disease cycles per gen -- tyring 2
# nb: in base case, no benefit of allele, but cost is fixed
# nb: disease ngens = 150

# so.... fig 1 will be comparision of structure across base model
cases <- expand.grid("transmission type" = 1:2, 
                     "compartments" = 1:3, # note SIR doesn't seem to have a big effect
                     "robustness" = 1:4) # 1 = mortality, 2 = transmission, 3 = recovery, 4 = demographic rescue only
# remove 1/3 compartment/robustness combo
cases <- cases[c(1:12, 15:24), ]
sim_dat <- NULL # this is where everything will get stored

# doing parallel
cl <- parallel::makeCluster(detectCores()-1)
doParallel::registerDoParallel(cl)

for (case in 19:22){
  
  print(case) # for error id and to check in
  
  case_vect <- base_vect
  # if transmission type is frequency (2), need to increase baseline transmission a bit
  if (cases$`transmission type`[case] == 2) {
    case_vect[6:8] <- 3 # transmission increased
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
    case_vect[18] <- 1
    case_vect[19] <- 0.5
  } # 1 = recovery, nb: recovery higher
  
  # 4 = no robustness, only change tranmission and such...
  # if (cases$robustness[case] == 4) {
  #   
  #   case_vect <- base_vect # no robustness
  #   
  #   # # demographics only...
  #   # case_vect[22:24] <- 2 # no cost?
  #   
  # } # make sure that everything is reset...
  
  # modify base parameters
  case_vect[3] <- cases$`transmission type`[case] # since numbers have meaning
  case_vect[1] <- cases$compartments[case] # diddo
  
  # run simulation
  case_dat <- run_gens(case_vect, 150, 500)
  
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

# save things!
saveRDS(fig1_dat, file = "figures/figure_data/fig1_dat.rds")

# colnames(sim_dat)[19] <- "init_R"
# get summaries to plot
prob_dat <- sim_dat %>% group_by(compartments, `transmission type`, robustness) %>% 
  summarize(`P(extinct)` = sum(extinct)/n(), #,
            `P(decline)` = sum(pop_drop50)/n()) # note for frequency dependent, force of infection is much much lower... think about how to make these more equalivalent
            # `P(population\nrecovery)` = sum(at_K95)/n(), # did the population recover?
            # `P(allele > 0.4)` = sum(r_allele_peak, na.rm = T)/n(), 
            # `P(lost disease)` = length(which(final_inf_prev == 0))/n())
er_dat <- sim_dat %>% filter(extinct == 0 & pop_drop50 == 1 & at_K95 == 1 & r_allele_peak45 == 1) %>% 
  group_by(compartments, `transmission type`, robustness) %>% 
  summarize(`P(rescue:\nevolutionary)` = n()/500)
dr_dat <- sim_dat %>% filter(extinct == 0 & pop_drop50 == 1 & at_K95 == 1 & r_allele_peak45 == 0) %>% 
  group_by(compartments, `transmission type`, robustness) %>% 
  summarize(`P(rescue:\ndemographic)` = n()/500)
# merge
tmp <- merge(prob_dat, er_dat, by = c("compartments", "transmission type", "robustness"), all = T)
plot_dat <- merge(tmp, dr_dat, by = c("compartments", "transmission type", "robustness"), all = T)
# rename
plot_dat$compartments <- recode(plot_dat$compartments, "1" = "SIX", "2" = "SIS", "3" = "SIR")
plot_dat$`transmission type` <- recode(plot_dat$`transmission type`, "1" = "density", "2" = "frequency")
plot_dat$robustness <- recode(plot_dat$robustness, "1" = "\u03bc", "2" = "\u03b2", "3" = "\u03d2", "4" = "N")
# convert to long -- note 4:10 for all the things
plot_long <- pivot_longer(plot_dat, cols = 4:7, names_to = "outcome", values_to = "value")
plot_long$value <- ifelse(is.na(plot_long$value), 0, plot_long$value)
# points plot -- also show varience? or issue with number of simulations?
fig1 <- ggplot(plot_long, aes(robustness, value, col = `transmission type`)) + 
  geom_hline(data = plot_long %>% filter(robustness == "N"), aes(yintercept = value), col = "gray70", lty = "dashed") + 
  geom_point(size = 3) + scale_color_manual(values = c("#ac1457", "#f1c4a2")) +
  facet_grid(rows = vars(outcome), cols = vars(compartments)) + 
  labs(x = "Evolutionary benefit", y = "Probability", col = "Transmission type") + 
  theme_bw()  + 
  theme(text = element_text(size = 16), legend.position = "bottom") 
# for box plots, need full data set!
tmp <- sim_dat %>% filter(extinct == 0 & pop_drop80 == 1 & at_K95 == 1 & r_allele_peak15 == 1 & robustness != 4) 
sim_long <- pivot_longer(tmp, cols = 1:18, names_to = "outcome", values_to = "value") %>% 
  filter(outcome %in% c("final_r_allele", # "max_r_allele", 
                        "final_inf_prev", # "max_inf_prev", 
                        # "final_pop_size", 
                        "firstK95")) # since only ER cases, don't really care about pop outcomes? should get time first at K
sim_long$compartments <- recode(sim_long$compartments, "1" = "SIX", "2" = "SIS", "3" = "SIR")
sim_long$`transmission type` <- recode(sim_long$`transmission type`, "1" = "density", "2" = "frequency")
sim_long$robustness <- recode(sim_long$robustness, "1" = "\u03bc", "2" = "\u03b2", "3" = "\u03d2", "4" = "N")
sim_long$outcome <- recode(sim_long$outcome, "final_r_allele" = "Final allele freq.", "max_r_allele" = "Max allele freq.", 
                                             "final_inf_prev" = "Final infection prev.", "max_inf_prev" = "Max infection prev.", 
                                             "final_pop_size" = "Final population size", "first_K95" = "Time to ER")
ggplot(sim_long %>% filter(compartments == "SIX"), aes(x = robustness, y = value, col = `transmission type`)) + 
  # geom_hline(data = plot_long %>% filter(robustness == "N"), aes(yintercept = median(value)), col = "gray70", lty = "dashed") + 
  geom_boxplot() + scale_color_manual(values = c("#ac1457", "#f1c4a2")) +
  facet_wrap(~outcome, scales = "free_y", nrow = 2) + 
  labs(x = "evolutionary benefit", y = NULL) + 
  theme_bw() + theme(text = element_text(size = 12), legend.position = "bottom") + 
  coord_cartesian(ylim = c(0, NA))
ggplot(sim_long %>% filter(compartments == "SIS"), aes(x = robustness, y = value, col = `transmission type`)) + 
  # geom_hline(data = plot_long %>% filter(robustness == "N"), aes(yintercept = median(value)), col = "gray70", lty = "dashed") + 
  geom_boxplot() + scale_color_manual(values = c("#ac1457", "#f1c4a2")) +
  facet_wrap(~outcome, scales = "free_y", nrow = 1) + 
  labs(x = "Evolutionary benefit", y = NULL, col = "Transmission type") + 
  theme_bw() + theme(text = element_text(size = 16), legend.position = "bottom") + 
  coord_cartesian(ylim = c(0, NA))

# comparing different cost/benefits
# so.... fig 1 will be comparision of structure across base model
sis_cb <- expand.grid("benefit" = c(0, 0.25, 0.5, 1, 1.25), # 0 = total benefit, 1.5 = some cost for B & M, but opposite for G
                      "cost" = c(1, 2), # 1 = some cost, 2 = no cost
                      "case" = c("\u03b2", "\u03bc", "\u03d2")) 
cb_dat <- NULL # this is where everything will get stored

# doing parallel
cl <- parallel::makeCluster(detectCores())
doParallel::registerDoParallel(cl)

for (cb in 11:20){
  
  print(cb) # for error id and to check in
  
  cb_vect <- base_vect
  cb_vect[1] <- 2 # SIS, with density
  # benefit update
  if (sis_cb$case[cb] == "\u03b2") {
    cb_vect[6] <- sis_cb$benefit[cb] # sis_cb$benefit[cb]*cb_vect[8]
    cb_vect[7] <- (cb_vect[6]+cb_vect[8])/2}
  if (sis_cb$case[cb] == "\u03bc") { 
    cb_vect[13] <- sis_cb$benefit[cb] # sis_cb$benefit[cb]*cb_vect[15]
    cb_vect[14] <- (cb_vect[13]+cb_vect[15])/2}
  if (sis_cb$case[cb] == "\u03d2") { 
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

fig2_dat <- cb_dat # G = 0.05

# save things!
saveRDS(fig2_dat, file = "figures/figure_data/fig2_dat.rds")

# get summaries to plot -- P(ER) + T_K + max/final R/inf (6) for the ER cases only
er_dat <- cb_dat %>% filter(extinct == 0 & pop_drop50 == 1 & at_K95 == 1 & r_allele_peak45 == 1) %>% 
  group_by(benefit, cost, case) %>% 
  summarize(`P(rescue: evolutionary)` = n()/500, 
            `final R allele` = mean(final_r_allele), 
            `max R allele` = mean(max_r_allele),
            `final inf prev` = mean(final_inf_prev),
            `max inf prev` = mean(max_inf_prev), 
            `time to K` = mean(firstK95))
cb_summ_long <- pivot_longer(er_dat, cols= 4:9, names_to = "outcome", values_to = "value")
cb_all_long <- pivot_longer(cb_dat %>% filter(extinct == 0 & pop_drop50 == 1 & at_K95 == 1 & r_allele_peak45 == 1) %>% select(c(8, 10, 18:21)), 
                            cols = 1:3, names_to = "outcome", values_to = "value")
ggplot(NULL, aes(x=as.factor(benefit), y=value, col = as.factor(cost))) + 
  geom_boxplot(data = cb_all_long) + 
  geom_point(data = cb_summ_long %>% filter(outcome == "P(rescue: evolutionary)"), size = 3) + scale_color_manual(values = c("#ac1457", "#f1c4a2")) + 
  facet_grid(rows = vars(outcome), cols = vars(case), scales = "free_y") + 
  theme_bw() + theme(text = element_text(size = 12), legend.position = "bottom") + 
  labs(x = "allele strength", col = "allele\nfecundity", y = NULL)

# then will do time series plot (note need a slightly different funciton to create these)

# other variable to check in the number of disease cycles between reproductions
sis_dc <- expand.grid("disease cycles" = 1:3, 
                      "robustness" = c(1, 3), # M, G only for now 
                      "event order" = c(2, 4, 6), # just for variation 
                      "transmission" = 2) # prelim suggests transmission types are similar
dc_dat <- NULL # this is where everything will get stored

# doing parallel
cl <- parallel::makeCluster(detectCores())
doParallel::registerDoParallel(cl)

for (case in 10:18){
  
  print(case) # for error id and to check in

  case_vect <- base_vect
  case_vect[1] <- 2 # SIS
  
  # if transmission type is frequency (2), need to increase baseline transmission a bit
  if (sis_dc$`transmission`[case] == 2) {
    case_vect[6:8] <- 2 # transmission increased
  }
  
  if (sis_dc$robustness[case] == 1) {
    case_vect[13] <- 0
    case_vect[14] <- 0.5
  } # 1 = mortality
  if (sis_dc$robustness[case] == 2) {
    case_vect[6] <- 0*case_vect[8] # since will be scaled by WW rate
    case_vect[7] <- 0.5*case_vect[8]
  } # 1 = transmission
  if (sis_dc$robustness[case] == 3) {
    case_vect[18] <- 1
    case_vect[19] <- 0.5
  } # 1 = recovery, nb: recovery higher
  
  # modify base parameters
  case_vect[3] <- sis_dc$`transmission`[case] # since numbers have meaning
  case_vect[2] <- sis_dc$`event order`[case] # diddo
  case_vect[29] <- sis_dc$`disease cycles`[case]
  
  # cost update -- not dependent on the case and always no cost
  case_vect[22] <- 2
  case_vect[23] <- (case_vect[22]+case_vect[24])/2
  
  # cycles update: equal time steps
  case_vect[4] <- floor(15/sis_dc$`disease cycles`[case])
  case_vect[5] <- floor(10/sis_dc$`disease cycles`[case])
  ngen_cycle <- floor(70/sis_dc$`disease cycles`[case])
  
  # run simulation
  case_dat <- run_gens(case_vect, ngen_cycle, 500)
  
  # save connect to case info
  case_dat$`disease cycles` <- rep(sis_dc$`disease cycles`[case])
  case_dat$`event order` <- rep(sis_dc$`event order`[case])
  case_dat$`transmission` <- rep(sis_dc$`transmission`[case])
  case_dat$`robustness` <- rep(sis_dc$`robustness`[case])
  
  # then append to sim_dat
  dc_dat <- rbind(dc_dat, case_dat)
  
}

# end parallel cluster
stopCluster(cl)

fig3_dat <- dc_dat # G = 0.05

# save things!
saveRDS(fig3_dat, file = "figure_data/fig3_dat.rds")

# get summaries to plot -- P(ER) + T_K + max/final R/inf (6) for the ER cases only --> note more generous def of pop drop
er_dat <- dc_dat %>% filter(extinct == 0 & pop_drop85 == 1 & at_K95 == 1 & r_allele_peak == 1) %>% 
  group_by(`disease cycles`, `robustness`, `event order`, `transmission`) %>% 
  summarize(`P(rescue:\nevolutionary)` = n()/500, 
            `final R allele` = mean(final_r_allele), 
            `max R allele` = mean(max_r_allele),
            `final inf prev` = mean(final_inf_prev),
            `max inf prev` = mean(max_inf_prev), 
            `time to K` = mean(first_K95))
dc_summ_long <- pivot_longer(er_dat, cols= 5:10, names_to = "outcome", values_to = "value")
dc_all_long <- pivot_longer(dc_dat %>% filter(extinct == 0 & pop_drop85 == 1 & at_K95 == 1 & r_allele_peak == 1) %>% select(c(6, 8, 16:19)), 
                             cols = 1:2, names_to = "outcome", values_to = "value")
ggplot(NULL, aes(x=as.factor(`disease cycles`), y=value, col=as.factor(robustness))) + 
  geom_boxplot(data = dc_all_long) + 
  geom_point(data = dc_summ_long %>% filter(outcome == "P(rescue: evolutionary)"), size = 3) + scale_color_manual(values = c("#ac1457", "#f1c4a2", "gray70")) + 
  facet_grid(row = vars(outcome), col = vars(`event order`), scales = "free_y") + 
  theme_bw() + theme(text = element_text(size = 12)) + labs(x = "Disease cycles between reproduction", col = NULL, y = NULL)
# think there's something weird with the rates, and if disease sticks around. might need to up the number of diseased individuals? 
