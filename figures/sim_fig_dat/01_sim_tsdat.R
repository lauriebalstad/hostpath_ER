library(tidyverse)
library(ggplot2)
library(viridis)

# so here really just want to iterate through a bunch of examples with different parameters
# grab random 10 parameters, no repeating?
# plus one that leads to coexistance
# plus one that has back selection... if possible

# 1159, 1133, 1172
# 1563, 1580, 1548, 1611

parm_mat <- readRDS("sensitivity/mat_var.rds")
grab_pars <- parm_mat[c(313:316, 1:3, 6, 1159, 1611), ] # NOTE: from older version of mat_var.rds...
coexist_pars <- c(1, 2, 2, # SIS, BGM, freq
                  100, 0.01, # timing things
                  0, 0, 3, 0.05, # 6-9: transmission
                  0.05, 0.05, 0.05, # background mort + 16
                  0, 1, 2, # 13-15 + 17: disease mort
                  0.01, 0.02, # mort sd
                  0.05, 0.05, 0.05, 0.01, # 18-21: recovery, unused bc SIX (1)
                  2, 2, 2, 0.2, 0.01, # reproduction & mutation
                  130, 4, # carrying capacity things 
                  1, 100, 9999)
back_pars <- c(2, 2, 1, # SIX, BGM, density
               100, 0.01, # timing things
               0, 0, 1, 0.05, # 6-9: transmission
               0.05, 0.05, 0.05, # background mort + 16
               1, 1, 1, # 13-15 + 17: disease mort
               0.01, 0.02, # mort sd
               0.05, 0.05, 0.05, 0.01, # 18-21: recovery, unused bc SIX (1)
               0.5, 1.25, 2, 0.2, 0.01, # reproduction & mutation
               130, 4, # carrying capacity things 
               1, 100, 9999) # nubmer of disease cycles per gen -- tyring 2
# combine
ts_pars <- rbind(grab_pars, coexist_pars, back_pars)
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

for (i in 1:dim(ts_pars)[1]) {
  
  ngens <- floor(150/ts_pars[i, 29])
  ts_df <- ts_data(ts_pars[i, ], ngens)
  ts_df$trial <- rep(i, dim(ts_df)[1])
  ts_df$ts <- ts_df$ts*ts_pars[i, 29] # rescale
  
  all_ts_dat <- rbind(all_ts_dat, ts_df)
  
}

all_ts_dat_tmp <- all_ts_dat # for saving

# get only the time steps around disease intro: ~ 100 -- 200
all_ts_dat <- all_ts_dat %>% filter(ts < 180 & ts > 90)
all_ts_dat$ts <- all_ts_dat$ts - 100 # 100 pre-disease intro steps
# add column for prevelence
all_ts_dat$prev <- all_ts_dat$I_pop/all_ts_dat$tot_pop
all_ts_dat$scaled_pop <- all_ts_dat$tot_pop/all_ts_dat$K_val
colnames(all_ts_dat) <- c("ts", "S_pop", "I_pop", "R_pop", "tot_pop", "K_val", "B - Adaptive allele\nfrequency", "trial", "C - Infection\nprevelence", "A - Scaled\npopulation size")

# save this time series run
saveRDS(all_ts_dat, file = "figures/figure_data/all_ts_dat.rds")
