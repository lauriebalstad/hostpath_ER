library(tidyverse)
library(ggplot2)
library(viridis)
source("funcs/00_model_timeseries.R")

set.seed(33) # repeatability

# so here really just want to iterate through a bunch of examples with different parameters
# note some combinations never result in population decline
# coexist_parsA <- c(2, 1, # SIS, BGM, freq
#                    0, 0, 0.7, 0.01, # 3-6 enviro dep transmission
#                    0, 0, 0.7, 0.05, # 7-10 dens dep transmission
#                    0.05, 0.05, 0.05, # 11-13 + 17 background mort 
#                    0, 1, 2, # 14-16 + 18 disease mort
#                    0.01, 0.02, # 17-18 mort sd
#                    0.5, 0.05, 0.05, 0.01, # 19-22 recovery, unused bc SIX (1)
#                    1, 1, 1, 0.001, # 23-27 reproduction & mutation
#                    130, 4, # 28-29 carrying capacity things 
#                    100, 1, 0.01, 0.1, 100) # 30-32 timing things and inital R, ngens placeholders
# coexist_parsB <- c(3, 4,
#                    -1, -1, -1, 0.01, # 3-6 enviro dep transmission
#                    0, 0, 1.5, 0.05, # 7-10 dens dep transmission
#                    0.05, 0.05, 0.05, # 11-13 + 17 background mort 
#                    0, 0.6, 1.2, # 14-16 + 18 disease mort
#                    0.01, 0.02, # 17-18 mort sd
#                    0.5, 0.05, 0.05, 0.01, # 19-22 recovery, unused bc SIX (1)
#                    0.5, 0.75, 1, 0.001, # 23-27 reproduction & mutation
#                    130, 4, # 28-29 carrying capacity things 
#                    100, 1, 0.01, 0.1, 100) # 30-32 timing things and inital R, ngens placeholders
# back_pars <- c(2, 1, 
#                -1, -1, -1, 0, 
#                0, 0, 1, 0.05, 
#                0.05, 0.05, 0.05, 
#                1.1, 1.1, 1.1, 
#                0.01, 0.02, 
#                0.05, 0.05, 0.05, 0.01, 
#                0.2, 0.2, 1, 0.001, 
#                130, 4, 
#                100, 1, 0.01, 0.1, 100) 
# ext_pars <- c(1, 2, 
#               -1, -1, -1, 0.1, 
#               0, 0, 1, 0.05, 
#               0.05, 0.05, 0.05, 
#               0, 1.6, 1.6, 
#               0.01, 0.02,
#               0.05, 0.05, 0.05, 0.01,
#               0.5, 0.75, 1, 0.001, 
#               130, 4, 
#               100, 1, 0.01, 0.1, 100) 
# demo_parsA <- c(3, 5, 
#                 -1, -1, -1, 0, 
#                 1, 1, 1, 0.05, # end: 10
#                 0.05, 0.05, 0.05, 
#                 0, 0, 1.2, 
#                 0.01, 0.02,
#                 0.5, 0.05, 0.05, 0.01, # start 19
#                 0.4, 1, 1, 0.005, 
#                 130, 4, 
#                 100, 1, 0.01, 0.1, 100) 
# demo_parsB <- c(1, 2, 
#                 -1, -1, -1, 0, 
#                 1, 1, 1, 0.05, # end: 10
#                 0.05, 0.05, 0.05, 
#                 0, 1.2, 1.2, 
#                 0.01, 0.02,
#                 1, 0.05, 0.05, 0.01, # start 19
#                 0.3, 0.65, 1, 0.005, 
#                 130, 4, 
#                 100, 1, 0.01, 0.1, 100) 

coexist_parsA <- c(2, 4, # compartments, event order
                   0.073/3, 0.073/3, 0.073/3, 0.001, # 3-6 no enviro dep transmission
                   0.073, 0.073, 0.073, 0.01, # 7-10 always dens dep transmission
                   0.0002, 0.0002, 0.0002, # 11-13 + 17 background mort 
                   0, 0.25, 0.5, # 14-16 + 18 disease mort
                   0.00005, 0.05, # 17-18 mort sd
                   0.005, 0.005, 0.005, 0.0005, # 19-22 recovery, unused bc SIX (1)
                   exp(0.0866)*1, exp(0.0866)*1, exp(0.0866), 0.0005, # 23-27 reproduction & mutation
                   # 0.6, 0.8, 1, 0.05, 0.005, # 23-27 reproduction & mutation
                   200, 4, # 28-29 carrying capacity things
                   100, 4, # 30-32 timing things # note change from d = 1 to d = 2???
                   0.01, 0.05, 80) # init R and init disease, ngens
coexist_parsB <- c(1, 4, # compartments, event order
                   0.073/3, 0.073/3, 0.073/3, 0.001, # 3-6 no enviro dep transmission
                   0.073, 0.073, 0.073, 0.01, # 7-10 always dens dep transmission
                   0.0002, 0.0002, 0.0002, # 11-13 + 17 background mort 
                   0, 0, 0.5, # 14-16 + 18 disease mort
                   0.00005, 0.05, # 17-18 mort sd
                   0.005*20, 0.005*20, 0.005, 0.0005, # 19-22 recovery, unused bc SIX (1)
                   exp(0.0866)*1, exp(0.0866)*1, exp(0.0866), 0.0005, # 23-27 reproduction & mutation
                   # 0.6, 0.8, 1, 0.05, 0.005, # 23-27 reproduction & mutation
                   200, 4, # 28-29 carrying capacity things
                   100, 4, # 30-32 timing things # note change from d = 1 to d = 2???
                   0.01, 0.05, 80) # init R and init disease, ngens
back_pars_A <- c(2, 1, # compartments, event order
                 0, 0.073/5/2, 0.073/5, 0.001, # 3-6 no enviro dep transmission
                 0, 0.073/2, 0.073, 0.01, # 7-10 always dens dep transmission
                 0.0002, 0.0002, 0.0002, # 11-13 + 17 background mort 
                 0, 0.25, 0.5, # 14-16 + 18 disease mort
                 0.00005, 0.05, # 17-18 mort sd
                 0.005*20, 0.005*10, 0.005, 0.0005, # 19-22 recovery, unused bc SIX (1)
                 exp(0.0866)*0.9, exp(0.0866)*0.95, exp(0.0866), 0.0005, # 23-27 reproduction & mutation
                 # 0.6, 0.8, 1, 0.05, 0.005, # 23-27 reproduction & mutation
                 200, 4, # 28-29 carrying capacity things
                 100, 4, # 30-32 timing things # note change from d = 1 to d = 2???
                 0.01, 0.05, 80) # init R and init disease, ngens
back_pars_B <- c(2, 6, # compartments, event order
                 -2, -2, -2, 0.001, # 3-6 no enviro dep transmission
                 0, 0, 0.073, 0.01, # 7-10 always dens dep transmission
                 0.0002, 0.0002, 0.0002, # 11-13 + 17 background mort 
                 0.4, 0.4, 0.4, # 14-16 + 18 disease mort
                 0.00005, 0.05, # 17-18 mort sd
                 0.005*40, 0.005*40, 0.005, 0.0005, # 19-22 recovery, unused bc SIX (1)
                 exp(0.0866)*0.3, exp(0.0866)*0.3, exp(0.0866), 0.0005, # 23-27 reproduction & mutation
                 # 0.6, 0.8, 1, 0.05, 0.005, # 23-27 reproduction & mutation
                 200, 4, # 28-29 carrying capacity things
                 100, 4, # 30-32 timing things # note change from d = 1 to d = 2???
                 0.01, 0.05, 80) # init R and init disease, ngens
ext_pars <- c(1, 1, # compartments, event order
              0.073/3, 0.073/3, 0.073/3, 0.001, # 3-6 no enviro dep transmission
              0.073*2, 0.073*2, 0.073*2, 0.01, # 7-10 always dens dep transmission
              0.0002*10, 0.0002*10, 0.0002*10, # 11-13 + 17 background mort 
              0, 0.4, 0.8, # 14-16 + 18 disease mort
              0.00005, 0.05, # 17-18 mort sd
              0.005, 0.005, 0.005, 0.0005, # 19-22 recovery, unused bc SIX (1)
              exp(0.0866)*0.6, exp(0.0866)*0.8, exp(0.0866), 0.0005, # 23-27 reproduction & mutation
              # 0.6, 0.8, 1, 0.05, 0.005, # 23-27 reproduction & mutation
              200, 4, # 28-29 carrying capacity things
              100, 4, # 30-32 timing things # note change from d = 1 to d = 2???
              0.01, 0.05, 80) # init R and init disease, ngens
demo_parsA <- c(2, 4, # compartments, event order
                -2, -2, -2, 0.001, # 3-6 no enviro dep transmission
                0, 0.073/2, 0.073, 0.01, # 7-10 always dens dep transmission
                0.0002, 0.0002, 0.0002, # 11-13 + 17 background mort 
                0.5, 0.5, 0.5, # 14-16 + 18 disease mort
                0.00005, 0.05, # 17-18 mort sd
                0.005, 0.005, 0.005, 0.0005, # 19-22 recovery, unused bc SIX (1)
                exp(0.0866)*0.4, exp(0.0866)*0.7, exp(0.0866), 0.0005, # 23-27 reproduction & mutation
                # 0.6, 0.8, 1, 0.05, 0.005, # 23-27 reproduction & mutation
                200, 4, # 28-29 carrying capacity things
                100, 4, # 30-32 timing things # note change from d = 1 to d = 2???
                0.01, 0.05, 80) # init R and init disease, ngens
demo_parsB <- c(3, 2, # compartments, event order
                -2, -2, -2, 0.001, # 3-6 no enviro dep transmission
                0, 0.073/2, 0.073, 0.01, # 7-10 always dens dep transmission
                0.0002, 0.0002, 0.0002, # 11-13 + 17 background mort 
                0.5, 0.5, 0.5, # 14-16 + 18 disease mort
                0.00005, 0.05, # 17-18 mort sd
                0.005, 0.005, 0.005, 0.0005, # 19-22 recovery, unused bc SIX (1)
                exp(0.0866)*1, exp(0.0866)*1, exp(0.0866), 0.0005, # 23-27 reproduction & mutation
                # 0.6, 0.8, 1, 0.05, 0.005, # 23-27 reproduction & mutation
                200, 4, # 28-29 carrying capacity things
                100, 4, # 30-32 timing things # note change from d = 1 to d = 2???
                0.01, 0.05, 80) # init R and init disease, ngens
no_dens_pars <- c(2, 4, # compartments, event order
                  Inf, Inf, Inf, 0.001, # 3-6 no enviro dep transmission
                  Inf, Inf, Inf, 0.01, # 7-10 always dens dep transmission
                  0.0002, 0.0002, 0.0002, # 11-13 + 17 background mort 
                  0, 0.25, 0.5, # 14-16 + 18 disease mort
                  0.00005, 0.05, # 17-18 mort sd
                  -2, -2, -2, 0.0005, # 19-22 recovery, unused bc SIX (1)
                  exp(0.0866)*9, exp(0.0866)*0.95, exp(0.0866), 0.0005, # 23-27 reproduction & mutation
                  # 0.6, 0.8, 1, 0.05, 0.005, # 23-27 reproduction & mutation
                  200, 4, # 28-29 carrying capacity things
                  100, 4, # 30-32 timing things # note change from d = 1 to d = 2???
                  0.01, 1, 60) # init R and init disease, ngens

# combine
ts_pars <- rbind(coexist_parsA, # coexist_parsA, 
                 coexist_parsB, 
                 ext_pars, ext_pars,
                 back_pars_A, back_pars_A,
                 back_pars_B, back_pars_B,
                 demo_parsA, demo_parsA, 
                 demo_parsB, demo_parsB, 
                 # coexist_parsB, 
                 back_pars_A, ext_pars, 
                 no_dens_pars, no_dens_pars)
# ts_pars[, 29] <- rep(1) # normalize disease cycles

all_ts_dat <- data.frame(
  ts = NULL,
  S_pop = NULL,
  I_pop = NULL,
  R_pop = NULL,
  tot_pop = NULL,
  K_val = NULL,
  r_freq = NULL,
  trial = NULL
)

# dim(ts_pars)[1] -- end after this to avoid abort issue
for (i in 1:dim(ts_pars)[1]) {
  
  # ts_pars[i, 31] <- 1 # normalize disease gens?
  # 1 disease gen per cycle...?
  ts_df <- ts_data(ts_pars[i, 1:33])
  ts_df$trial <- rep(i, dim(ts_df)[1])
  ts_df$ts <- 1:length(ts_df$ts) # rescale
  
  # if no pop drop, remove simulation, using generous cut off
  if (all(ts_df$tot_pop > ts_df$K_val*0.85)) ts_df <- NULL
  
  all_ts_dat <- rbind(all_ts_dat, ts_df)
  
}

all_ts_dat_tmp <- all_ts_dat # for saving

# get only the time steps around disease intro: ~ 100 -- 200
# all_ts_dat <- all_ts_dat %>% filter(ts > 90 & ts < 130) # ts < 300 & 
all_ts_dat$ts <- all_ts_dat$ts - 25 # 100 pre-disease intro steps
# add column for prevelence
all_ts_dat$prev <- ifelse(is.nan(all_ts_dat$I_pop/all_ts_dat$tot_pop), 0, all_ts_dat$I_pop/all_ts_dat$tot_pop)
# rescale pop to maximum of its size -- NOTE: different than dividing by K at that timestep
all_ts_dat <- all_ts_dat %>% group_by(trial) %>% mutate(max_pop = (max(tot_pop)))
all_ts_dat$scaled_pop <- all_ts_dat$tot_pop/all_ts_dat$max_pop
# remove NaN -- just 0s
all_ts_dat$r_freq <- ifelse(is.nan(all_ts_dat$r_freq), 0, all_ts_dat$r_freq)
colnames(all_ts_dat) <- c("ts", "S_pop", "I_pop", "R_pop", "tot_pop", "K_val", "Adaptive allele\nfrequency", "trial", "Infection\nprevelence", "max_pop", "Scaled\npopulation size")
first_zero <- all_ts_dat %>% group_by(trial) %>% summarize (first_zero = ts[which(tot_pop == 0)[1]]) # sometimes simulations add a few zeros after extinction--remove that
first_zero$first_zero <- ifelse(is.na(first_zero$first_zero), 80, first_zero$first_zero)
all_ts_dat <- merge(all_ts_dat, first_zero) %>% filter(ts <= first_zero)

# save this time series run
saveRDS(all_ts_dat, file = "dat/all_ts_dat_0602.rds")
