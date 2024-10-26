library(tidyverse)

# so here really just want to iterate through a bunch of examples with different parameters
# grab random 10 parameters, no repeating?
# plus one that leads to coexistance
# plus one that has back selection... if possible

# 1159, 1133, 1172
# 1563, 1580, 1548, 1611

grab_pars <- parm_mat[c(359:369), ] # 359:369 works
coexist_pars <- parm_mat[c(1159, 1611), ] # 1159, 1611
back_pars <- c(2, 2, 1, # SIX, BGM, density
               10, 10, # timing things
               0, 0, 3, 0.05, # 6-9: transmission
               0.05, 0.05, 0.05, # background mort + 16
               1, 1, 1, # 13-15 + 17: disease mort
               0.01, 0.02, # mort sd
               0.05, 0.05, 0.05, 0.01, # 18-21: recovery, unused bc SIX (1)
               1, 1, 2, 0.2, 0.01, # reproduction & mutation
               130, 4, # carrying capacity things 
               1) # nubmer of disease cycles per gen -- tyring 2
# combine
ts_pars <- rbind(grab_pars, coexist_pars, back_pars, back_pars, back_pars)

all_ts_dat <- data.frame(
  ts = NULL,
  S_pop = NULL,
  I_pop = NULL,
  R_pop = NULL,
  tot_pop = NULL,
  K_val = NULL,
  r_freq = NULL
)

for (i in 1:dim(ts_pars)[1]) {
  
  ngens <- floor(100/ts_pars[i, 29])
  ts_df <- ts_data(ts_pars[i, ], ngens)
  ts_df$trial <- rep(i, dim(ts_df)[1])
  
  all_ts_dat <- rbind(all_ts_dat, ts_df)
  
}

# add column for prevelence
all_ts_dat$prev <- all_ts_dat$I_pop/all_ts_dat$tot_pop

# save this time series run
# saveRDS(all_ts_dat, file = "figure_data/all_ts_dat.rds")

# rotate data
all_ts_long <- pivot_longer(all_ts_dat %>% select(c(1, 5, 7:9)), cols = c(2, 3, 5), names_to = "outcome", values_to = "value")
# grab back selection case -- 15
# grab extinction case -- 10
# grab coexistance case -- 13

ggplot(all_ts_long %>% filter(trial %in% c(4, 7:16)), aes(ts, value)) + # 7, 8, 9, 10 gives exinction
  geom_line(aes(group = trial), col = "gray80", lwd = 0.75) +
  geom_line(data = all_ts_long %>% filter(trial == 10), col = "#DB6341", lwd = 0.85) + # extinct
  geom_line(data = all_ts_long %>% filter(trial == 12), col = "#f1c4a2", lwd = 0.85) + # DR only, lost disease
  geom_line(data = all_ts_long %>% filter(trial == 4), col = "black", lwd = 0.85) + # coexistance
  geom_line(data = all_ts_long %>% filter(trial == 14), col = "#ac1457", lwd = 0.85) + # back select
  facet_wrap(~outcome, scale = "free_y") + theme_bw() + 
  labs(x = "generation", y = NULL)

