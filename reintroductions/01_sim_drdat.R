library(tidyverse)

# so first, need to simulate out a population with disease drop
coexist_pars <- c(2, 2, 0.1, # SIS, BGM, freq
                  100, 0.01, # timing things -- d0 and ngens; hoping to stop ngens at drop?
                  1, 1, 1, 0.05, # 6-9: transmission
                  0.05, 0.05, 0.05, # background mort + 16
                  1, 1, 1.2, # 13-15 + 17: disease mort
                  0.01, 0.02, # mort sd
                  1, 0.5, 0.05, 0.01, # 18-21: recovery, unused bc SIX (1)
                  0.5, 1.25, 2, 0.2, 0.01, # reproduction & mutation
                  130, 4, # carrying capacity things 
                  1, 100)

classic_pars <- c(2, 2, 0.1, # SIS, BGM, freq
                  100, 0.01, # timing things -- d0 and ngens; hoping to stop ngens at drop?
                  1.5, 2, 2, 0.05, # 6-9: transmission
                  0.05, 0.05, 0.05, # background mort + 16
                  0.4, 0.8, 1.2, # 13-15 + 17: disease mort
                  0.01, 0.02, # mort sd
                  0.25, 0.25, 0.05, 0.01, # 18-21: recovery, unused bc SIX (1)
                  2, 2, 2, 0.2, 0.01, # reproduction & mutation
                  130, 4, # carrying capacity things 
                  1, 100)

# run to get pre-reintro conditions
pre_ri <- dr_init(coexist_pars, 10) 
dim(pre_ri)

# then run dr_data
reintro_pars <- coexist_pars
reintro_pars[4] <- 0 # recall parameter 4's meaning is new
reintro_pars[5] <- 70 # recall parameter 5's meaning is new
reintro_pars[26] <- 0.01
# first get no reintro example
no_ri <- dr_data(reintro_pars, pre_ri, 0, 0)
large_ri <- dr_data(reintro_pars, pre_ri, 52, 1)
medium_ri  <- dr_data(reintro_pars, pre_ri, 26, 1)
small_ri  <- dr_data(reintro_pars, pre_ri, 52, 4)

# plot things
ymax <- max(c(no_ri$tot_pop, large_ri$tot_pop, medium_ri$tot_pop, small_ri$tot_pop))
par(mfrow = c(2, 2))
plot(no_ri$ts, no_ri$tot_pop, ylim = c(0, ymax))
lines(no_ri$ts, ymax*no_ri$r_freq, col = "blue"); lines(no_ri$ts, no_ri$I_pop, col = "red")
plot(large_ri$ts, large_ri$tot_pop, ylim = c(0, ymax))
lines(large_ri$ts, ymax*large_ri$r_freq, col = "blue"); lines(large_ri$ts, large_ri$I_pop, col = "red")
plot(medium_ri$ts, medium_ri$tot_pop, ylim = c(0, ymax))
lines(medium_ri$ts, ymax*medium_ri$r_freq, col = "blue"); lines(medium_ri$ts, medium_ri$I_pop, col = "red")
plot(small_ri$ts, small_ri$tot_pop, ylim = c(0, ymax))
lines(small_ri$ts, ymax*small_ri$r_freq, col = "blue"); lines(small_ri$ts, small_ri$I_pop, col = "red")

# iterate through a bunch of the cases and summarise
summ_dat_backsel <- NULL

for (i in c(1:400)) {
  
  no_ri <- dr_data(reintro_pars, pre_ri, 0, 0); no_ri$at_K <- (no_ri$tot_pop >= no_ri$K_val*0.95)
  large_ri <- dr_data(reintro_pars, pre_ri, 52, 1); large_ri$at_K <- (large_ri$tot_pop >= large_ri$K_val*0.95)
  medium_ri <- dr_data(reintro_pars, pre_ri, 26, 1); medium_ri$at_K <- (medium_ri$tot_pop >= medium_ri$K_val*0.95)
  small_ri  <- dr_data(reintro_pars, pre_ri, 52, 4); small_ri$at_K <- (small_ri$tot_pop >= small_ri$K_val*0.95)
  
  no_ri_summ <- data.frame(ex = ifelse(any(no_ri$tot_pop == 0), 1, 0), 
                           er = ifelse(any(tail(no_ri$at_K, 5)) & 
                                          max(no_ri$r_freq) > 0.4 & 
                                          all(no_ri$tot_pop > 0), 1, 0), 
                           dr = ifelse(any(tail(no_ri$at_K, 5)) & 
                                          max(no_ri$r_freq) < 0.4 & 
                                         all(no_ri$tot_pop > 0), 1, 0), 
                           time_er = which(no_ri$at_K)[1],
                           max_inf = max(no_ri$I_pop/no_ri$tot_pop, na.rm = T, nan.rm = T), 
                           final_inf = mean(tail(no_ri$I_pop/no_ri$tot_pop, 15), na.rm = T, nan.rm = T), 
                           final_r = mean(tail(no_ri$r_freq, 5), na.rm = T, nan.rm = T), 
                           final_pop = mean(tail(no_ri$tot_pop, 5), na.rm = T, nan.rm = T), 
                           strat = "no ri", case = i)
  lg_ri_summ <- data.frame(ex = ifelse(any(large_ri$tot_pop == 0), 1, 0), 
                           er = ifelse(any(tail(large_ri$at_K, 15)) & 
                                         max(large_ri$r_freq) > 0.4 & 
                                         all(large_ri$tot_pop > 0), 1, 0), 
                           dr = ifelse(any(tail(large_ri$at_K, 5)) & 
                                         max(large_ri$r_freq) < 0.4 & 
                                         all(large_ri$tot_pop > 0), 1, 0), 
                           time_er = which(large_ri$at_K)[1],
                           max_inf = max(large_ri$I_pop/large_ri$tot_pop, na.rm = T, nan.rm = T), 
                           final_inf = mean(tail(large_ri$I_pop/large_ri$tot_pop, 5), na.rm = T, nan.rm = T), 
                           final_r = mean(tail(large_ri$r_freq, 5), na.rm = T, nan.rm = T), 
                           final_pop = mean(tail(large_ri$tot_pop, 5), na.rm = T, nan.rm = T), 
                           strat = "large ri", case = i)
  md_ri_summ <- data.frame(ex = ifelse(any(medium_ri$tot_pop == 0), 1, 0), 
                           er = ifelse(any(tail(medium_ri$at_K, 5)) & 
                                         max(medium_ri$r_freq) > 0.4 & 
                                         all(medium_ri$tot_pop > 0), 1, 0), 
                           dr = ifelse(any(tail(medium_ri$at_K, 5)) & 
                                         max(medium_ri$r_freq) < 0.4 & 
                                         all(medium_ri$tot_pop > 0), 1, 0), 
                           time_er = which(medium_ri$at_K)[1],
                           max_inf = max(medium_ri$I_pop/medium_ri$tot_pop, na.rm = T, nan.rm = T), 
                           final_inf = mean(tail(medium_ri$I_pop/medium_ri$tot_pop, 5), na.rm = T, nan.rm = T), 
                           final_r = mean(tail(medium_ri$r_freq, 5), na.rm = T, nan.rm = T), 
                           final_pop = mean(tail(medium_ri$tot_pop, 5), na.rm = T, nan.rm = T), 
                           strat = "medium ri", case = i)
  sm_ri_summ <- data.frame(ex = ifelse(any(small_ri$tot_pop == 0), 1, 0), 
                           er = ifelse(any(tail(small_ri$at_K, 5)) & 
                                         max(small_ri$r_freq) > 0.4 & 
                                         all(small_ri$tot_pop > 0), 1, 0), 
                           dr = ifelse(any(tail(small_ri$at_K, 5)) & 
                                         max(small_ri$r_freq) < 0.4 & 
                                         all(small_ri$tot_pop > 0), 1, 0), 
                           time_er = which(small_ri$at_K)[1],
                           max_inf = max(small_ri$I_pop/small_ri$tot_pop, na.rm = T, nan.rm = T), 
                           final_inf = mean(tail(small_ri$I_pop/small_ri$tot_pop, 5), na.rm = T, nan.rm = T), 
                           final_r = mean(tail(small_ri$r_freq, 5), na.rm = T, nan.rm = T), 
                           final_pop = mean(tail(small_ri$tot_pop, 5), na.rm = T, nan.rm = T), 
                           strat = "small ri", case = i)
  
  summ_dat_backsel <- rbind(summ_dat_backsel, no_ri_summ, lg_ri_summ, md_ri_summ, sm_ri_summ)
  
}

backselection_DR <- list(coexist_pars, reintro_pars, pre_ri, summ_dat_backsel)
saveRDS(backselection_DR, file = "reintroductions/backselection_DR.rds")

# classic_DR <- list(classic_pars, reintro_pars, pre_ri, summ_dat_classic)
# saveRDS(classic_DR, file = "reintroductions/classic_DR.rds")

summ_dat_backsel$sit <- rep("high cost")
summ_dat_classic$sit <- rep("classic")
summ_dat_tot <- rbind(summ_dat_backsel, summ_dat_classic)

# get summary stats
summ_stat <- summ_dat_tot %>% group_by(strat, sit) %>%
  summarise(`P(ex)` = sum(ex)/length(ex), 
            `P(ER)` = sum(er)/length(er)) # dr too rare bc always some enviro transmission
er_dat <- summ_dat_tot %>% filter(er == 1)
summ_long <- pivot_longer(summ_stat, cols = 3:4)
all_long <- pivot_longer(er_dat %>% select(c(4, 7:11)), cols = 1:3)

ggplot(data = NULL, aes(strat, value, col = sit)) + 
  geom_point(data = summ_long, size = 3) + geom_boxplot(data = all_long) + 
  facet_wrap(~name, ncol = 1, scale = "free_y") + 
  scale_color_manual(values = c("#ac1457", "#f1c4a2")) +
  theme_bw() + labs(xlab = "reintroduction strategy", ylab = NULL, col = "Parameter case")
