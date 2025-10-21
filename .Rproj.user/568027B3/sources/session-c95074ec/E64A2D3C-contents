library(ggplot2)
library(cowplot)

# want to check how run_gens looks over different numbers of internal simulations
# might need to look also at how number of generations starts to matter
# get a sense of how many reps are needed

# starting with parameters that i know will kind of work....
test_vect <- c(1, 2, 1, # SIX, BGM, density
               15, 10, # timing things
               1, 1, 1, 0.05, # transmission
               0.05, 0.05, 0.05, # background mort
               0.1, 0.5, 1, # disease mort
               0.01, 0.02, # mort sd
               0.1, 0.1, 0.1, 0.05, # recovery, unused bc SIX (1) 
               1.5, 1.75, 2, 0.2, 0.01, # reproduction & mutation
               100, 4) # carrying capacity things 

# comp_str <- parm_vect[1]; event_order <- parm_vect[2]; trans_type <- parm_vect[3] # structural
# d_0 <- parm_vect[4]; r_0 <- parm_vect[5] # disease intro and allele intro
# b_RR <- parm_vect[6]; b_WR <- parm_vect[7]; b_WW <- parm_vect[8]; b_sd <- parm_vect[9] # transmission
# m_SRR <- parm_vect[10]; m_SWR <- parm_vect[11]; m_SWW <- parm_vect[12]; m_IRR <- parm_vect[13]; m_IWR <- parm_vect[14]; m_IWW <- parm_vect[15]; m_Ssd <- parm_vect[16]; m_Isd <- parm_vect[17] # mortality
# r_RR <- parm_vect[18]; r_WR <- parm_vect[19]; r_WW <- parm_vect[20]; r_sd <- parm_vect[21] # recovery
# l_RR <- parm_vect[22]; l_WR <- parm_vect[23]; l_WW  <- parm_vect[24]; l_sd <- parm_vect[25]; mut_rate <- parm_vect[26] # reproduction
# K <- parm_vect[27]; K_sd <- parm_vect[28] # carrying capacity
# n <- parm_vect[27] # number of individuals at start of simulation, all WW -- start at K

# doing parallel
cl <- parallel::makeCluster(detectCores())
doParallel::registerDoParallel(cl)

# ---- 30 gens ----
# running 500 sims 5 times
# make huge data frame to draw from
# use mat_var[71, ] --> does have extinction sometimes ~10/500?
# might need to modify mat_var, especially initial starting pop numbers to increase extinction risk?
# using modified version of model where output is the data frame
# start by testing 30 generations
r1 <- run_gens(test_vect, 70, 500)
r2 <- run_gens(test_vect, 70, 500)
r3 <- run_gens(test_vect, 70, 500)
r4 <- run_gens(test_vect, 70, 500)
r5 <- run_gens(test_vect, 70, 500)
r6 <- run_gens(test_vect, 70, 500)
r7 <- run_gens(test_vect, 70, 500)
r8 <- run_gens(test_vect, 70, 500)
r9 <- run_gens(test_vect, 70, 500)
r10 <- run_gens(test_vect, 70, 500)
r11 <- run_gens(test_vect, 70, 1000)
r12 <- run_gens(test_vect, 70, 1000)
r13 <- run_gens(test_vect, 70, 1000)
r14 <- run_gens(test_vect, 70, 1000)
r15 <- run_gens(test_vect, 70, 1000)
r16 <- run_gens(test_vect, 70, 1000)
r_all_70 <- rbind(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, 
                  r11, r12, r13, r14, r15) 
r_all_70$num <- 1:dim(r_all_70)[1]

saveRDS(r_all_70, file = "sensitivity_data/r_all_70.rds")

# end parallel cluster
stopCluster(cl)


# overall: now want to get a sense of...
# P(ex)
sum(r_all_70$extinct)/length(r_all_70$extinct) # 0.4654
p_ext <- ggplot(data = r_all_70, aes(x = extinct)) + 
  geom_density(col = "#FF61CC") + geom_vline(aes(xintercept = sum(r_all_70$extinct)/length(r_all_70$extinct)),lty = "dashed") +
  theme_classic()
# P(ER -- pop_drop + at_K + allele_peak) 
length(which(r_all_70$extinct == 0 &
               r_all_70$pop_drop == 1 & 
               r_all_70$at_K == 1 & 
               r_all_70$r_allele_peak == 1))/dim(r_all_70)[1] # 0.4998
er_dat <- r_all_70 %>% filter(extinct == 0 & 
                                pop_drop == 1 &
                                at_K == 1 &
                                r_allele_peak == 1)
prob_df <- data.frame(num = 1:dim(r_all_70)[1],
  p_er = ifelse((r_all_70$extinct == 0 &
                              r_all_70$pop_drop == 1 & 
                              r_all_70$at_K == 1 & 
                              r_all_70$r_allele_peak == 1), 1, 0), 
  p_atK = ifelse((r_all_70$extinct == 0 & r_all_70$at_K), 1, 0))
p_er <- ggplot(data = prob_df, aes(x = p_er)) + 
  geom_density(col = "#C77CFF") + 
  geom_vline(aes(xintercept = length(which(r_all_70$extinct == 0 &
                                             r_all_70$pop_drop == 1 & 
                                             r_all_70$at_K == 1 & 
                                             r_all_70$r_allele_peak == 1))/dim(r_all_70)[1]), lty = "dashed") + 
  theme_classic()
p_atK <- ggplot(data = prob_df, aes(x = p_atK)) + 
  geom_density(col = "#00A9FF") + 
  geom_vline(aes(xintercept = length(which(r_all_70$extinct == 0 &
                                             r_all_70$at_K == 1))/dim(r_all_70)[1]), lty = "dashed") + 
  theme_classic()
# --> follow up with avg infection at end
mean(er_dat$final_inf_prev) # 0.22... def bimodal tho
er_finalINF <- ggplot(data = er_dat, aes(x = final_inf_prev)) + 
  geom_density(col = "#F8766D") + geom_vline(aes(xintercept = mean(final_inf_prev)), lty = "dashed") + 
  theme_classic() 
er_maxINF <- ggplot(data = er_dat, aes(x = max_inf_prev)) + 
  geom_density(col = "#7CAE00") + geom_vline(aes(xintercept = mean(max_inf_prev)), lty = "dashed") + 
  theme_classic()
# --> follow up with time to reach K/allele_peak
mean(er_dat$first_K) # 43.45 ish... pretty centered, some long tail
er_timeK <- ggplot(data = er_dat, aes(x = first_K)) + 
  geom_density(col = "#00BFC4") + geom_vline(aes(xintercept = mean(first_K)), lty = "dashed") + 
  theme_classic() 
# --> follow up with final r allele & max r allele
mean(er_dat$max_r_allele) # 0.93... is bimodal. could make 0.4 cut point???
er_maxR <- ggplot(data = er_dat, aes(x = max_r_allele)) + 
  geom_density(col = "#00BE67") + geom_vline(aes(xintercept = mean(max_r_allele)), lty = "dashed") + 
  theme_classic() 
mean(er_dat$final_r_allele) # 0.28... def bimodal tho
er_finalR <- ggplot(data = er_dat, aes(x = final_r_allele)) + 
  geom_density(col = "#CD9600") + geom_vline(aes(xintercept = mean(final_r_allele)), lty = "dashed") + 
  theme_classic() 
# --> follow up wiht some comparision plots: 
plot(er_dat$final_r_allele, er_dat$final_inf_prev) # mostly positive but some mix... e.g, high r w/high prev and w/low prev
plot(er_dat$max_r_allele, er_dat$final_inf_prev) # very positive
# P(DR only -- pop_drop + at_K + NO allele_peak)
length(which(r_all_70$extinct == 0 &
               r_all_70$pop_drop == 1 & 
               r_all_70$at_K == 1 & 
               r_all_70$r_allele_peak == 0))/dim(r_all_70)[1] # 0.0336
dr_dat <- r_all_70 %>% filter(extinct == 0 & 
                                pop_drop == 1 &
                                at_K == 1 &
                                r_allele_peak == 0)
# --> follow up with r allele info
mean(dr_dat$max_r_allele) # 0.10... is bimodal 
ggplot(data = dr_dat, aes(x = max_r_allele)) + geom_density() + theme_classic() 
mean(dr_dat$final_r_allele) # 0.015... def bimodal tho
ggplot(data = dr_dat, aes(x = final_r_allele)) + geom_density() + theme_classic() 
# --> and disease info (i.e., did disease get lost?)
mean(dr_dat$final_inf_prev) # 0
mean(dr_dat$max_inf_prev) # 0.51. generally a peak but lots of heterogenity
ggplot(data = dr_dat, aes(x = max_inf_prev)) + geom_density() + theme_classic() 
# P(lost infection) 
length(which(r_all_70$final_inf_prev == 0 &
         r_all_70$extinct == 0))/dim(r_all_70)[1] # 0.40... pretty high
ni_dat <- r_all_70 %>% filter(extinct == 0, 
                              final_inf_prev == 0)
# --> follow up with avg allele + peak allele
# --> follow up with r allele info
mean(ni_dat$max_r_allele) # 0.84 does have good peak
ggplot(data = ni_dat, aes(x = max_r_allele)) + geom_density() + theme_classic() 
mean(ni_dat$final_r_allele) # 0.03, strong peak
# --> and disease info (i.e., did disease get lost?)
mean(ni_dat$max_inf_prev) # 0.94 -- pretty consistant with long tail
ggplot(data = ni_dat, aes(x = max_inf_prev)) + geom_density() + theme_classic() 

# now need to do the draws and such
num_draws <- c(10, 25, 50, 75, 100, 150, 200, 250, 300, 350, 400, 500, 750, 1000, 1500)
case_num <- 1:5
dat_plot <- expand.grid(num_draws = num_draws, case_num = case_num)
for (i in 1:dim(dat_plot)[1]) {
  
  # get row draws
  row_samp <- sample(r_all_70$num, dat_plot$num_draws[i], replace = F)
  # get filtered data frame
  r_tmp <- r_all_70 %>% filter(num %in% row_samp)
  # get key summary stats overall
  dat_plot$p_ext[i] <- sum(r_tmp$extinct)/dim(r_tmp)[1] # extinction
  dat_plot$p_er[i] <- length(which(r_tmp$extinct == 0 &
                                     r_tmp$pop_drop == 1 & 
                                     r_tmp$at_K == 1 & 
                                     r_tmp$r_allele_peak == 1))/dim(r_tmp)[1] # evolutionary rescue
  dat_plot$p_atK[i] <- length(which(r_tmp$at_K == 1 & r_tmp$extinct == 0))/dim(r_tmp)[1]
  # filter to only yes er
  er_tmp <- r_tmp %>% filter(extinct == 0 & pop_drop == 1 & at_K == 1 & r_allele_peak == 1)
  # get er stats
  dat_plot$er_maxR[i] <- mean(er_tmp$max_r_allele)
  dat_plot$er_finalR[i] <- mean(er_tmp$final_r_allele)
  dat_plot$er_maxINF[i] <- mean(er_tmp$max_inf_prev)
  dat_plot$er_finalINF[i] <- mean(er_tmp$final_inf_prev)
  dat_plot$er_timeK[i] <- mean(er_tmp$first_K)
  
}

# full sim values
# filter to only yes er
er_all <- r_all_70 %>% filter(extinct == 0 & pop_drop == 1 & at_K == 1 & r_allele_peak == 1)
full_dat <- data.frame(p_ext = sum(r_all_70$extinct)/dim(r_all_70)[1], # extinction
                       p_er = length(which(r_all_70$extinct == 0 &
                                            r_all_70$pop_drop == 1 & 
                                            r_all_70$at_K == 1 & 
                                            r_all_70$r_allele_peak == 1))/dim(r_all_70)[1], # evolutionary rescue
                       p_atK = length(which(r_all_70$at_K == 1 & r_all_70$extinct == 0))/dim(r_all_70)[1],
                       er_maxR = mean(er_all$max_r_allele), 
                       er_finalR = mean(er_all$final_r_allele), 
                       er_maxINF = mean(er_all$max_inf_prev),
                       er_finalINF = mean(er_all$final_inf_prev),
                       er_timeK = mean(er_all$first_K)
                       )

graph_dat <- reshape(dat_plot, 
                     varying = 3:10, v.names = "val", 
                     times = colnames(dat_plot)[3:10], timevar = "output",
                     direction = "long")

graph_dat_full <- reshape(full_dat, 
                     varying = 1:8, v.names = "val", 
                     times = colnames(full_dat)[1:8], timevar = "output",
                     direction = "long")

r70 <- ggplot(graph_dat, aes(num_draws, val, col = output)) + 
  geom_point(show.legend = F) + 
  facet_wrap(~output, scales = "free_y") + 
  theme_classic() + 
  geom_hline(data = graph_dat_full, aes(yintercept = val), lty = 2)

# r70

# and a cowplot of all the histograms to get a sense of normal-ness
all_dat_hists <- plot_grid(er_finalINF, er_finalR, er_maxINF, 
                           er_maxR, er_timeK, p_atK, 
                           p_er, p_ext, nrow = 3)
# all_dat_hists

plot_grid(r70, all_dat_hists, nrow = 1)

ggplot(r_all_70, aes(max_r_allele, 
                     final_inf_prev, 
                     col = final_r_allele,
                     shape = as.factor(extinct))) + 
  geom_point() + scale_color_gradient(low = "#00A9FF", high = "#C77CFF") + 
  geom_point(data = er_dat, pch = 1, col = "black") +
  theme_classic()

ggplot(r_all_70, aes(max_r_allele, 
                     max_inf_prev, 
                     col = final_r_allele,
                     shape = as.factor(extinct))) + 
  geom_point() + scale_color_gradient(low = "#00A9FF", high = "#C77CFF") + 
  geom_point(data = er_dat, pch = 1, col = "black") +
  theme_classic()

ggplot(r_all_70, aes(max_inf_prev, 
                     final_r_allele, 
                     col = final_inf_prev,
                     shape = as.factor(extinct))) + 
  geom_point() + scale_color_gradient(low = "#00A9FF", high = "#C77CFF") + 
  geom_point(data = er_dat, pch = 1, col = "black") +
  theme_classic()