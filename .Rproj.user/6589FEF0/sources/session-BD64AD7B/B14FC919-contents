library(dplyr)

# output order: 
# [1] p_extinct, p_demo_recovery, avg_time_K, sd_time_K,
# [5] avg_pop_size_tf, sd_pop_size_tf, avg_pop_size_given_NE, sd_pop_size_given_NE,
# [9] p_lost_I, avg_time_lost_I, sd_time_lost_I,
# [12] avg_infect_class_I_tf, sd_infect_class_I_tf, avg_infect_class_I_given_NE, sd_infect_class_I_given_NE,
# [16] avg_freq_R_tf, sd_freq_R_tf, avg_max_R, avg_time_max_R

# want to check how run_gens looks over different numbers of internal simulations
# might need to look also at how number of generations starts to matter
num_gens <- c(15, 30, 45) 
num_reps <- c(2, 10, 50, 100, 250, 400, 750)
test_gen_reps_dat <- expand.grid(gens=num_gens, reps=num_reps)
# general set up for other parameters
tmp_dat <- NULL
for (i in 1:10) {
tmp <- mapply(run_gens, 1, 2,
                        2, 0, 10, 0.2, 
                        1, 1, 1, 0.1, 
                        0.1, 0.1, 0.1, 0.1, 1, 1, 0.1, 0.1,
                        1, 1, 1, 0.1, 
                        2, 4, 4, 0.5, 
                        150, 20, 
              test_gen_reps_dat$gens, test_gen_reps_dat$reps)
# SIX structure with "B","G","M" event order
# complete mortality blocking, fitness cost
# note infection cost is only mortality; no fitness cost to infection
tmp_tmp <- cbind(test_gen_reps_dat, t(tmp))
tmp_dat <- rbind(tmp_dat, tmp_tmp)
}

tgr <- tmp_dat
colnames(tgr) <- c("gens", "reps", 
                   "p_extinct", "p_demo_recovery", "avg_time_K", "sd_time_K",
                   "avg_pop_size_tf", "sd_pop_size_tf", "avg_pop_size_given_NE", "sd_pop_size_given_NE",
                   "p_lost_I", "avg_time_lost_I", "sd_time_lost_I",
                   "avg_infect_class_I_tf", "sd_infect_class_I_tf", "avg_infect_class_I_given_NE", "sd_infect_class_I_given_NE",
                   "avg_freq_R_tf", "sd_freq_R_tf", "avg_max_R", "avg_time_max_R")
plot(tgr$reps, tgr$p_extinct)
plot(tgr$reps, tgr$avg_time_K)
plot(tgr$reps, tgr$avg_freq_R_tf)

extinct_plot <- ggplot(tgr, aes(unlist(reps), unlist(p_extinct), col = unlist(gens))) + 
  geom_point() + labs(x = "number of repetitions (bnum)", 
                      y = "P(extinct)", 
                      col = "number of \ngenerations \n(ngens)") + theme_classic() 
timeK_plot <- ggplot(tgr, aes(unlist(reps), unlist(avg_time_K), col = unlist(gens))) + 
  geom_point() + labs(x = "number of repetitions (bnum)", 
                      y = "time to reach K", 
                      col = "number of \ngenerations \n(ngens)") + theme_classic() 
lostI_plot <- ggplot(tgr, aes(unlist(reps), unlist(p_lost_I), col = unlist(gens))) + 
  geom_point() + labs(x = "number of repetitions (bnum)", 
                      y = "P(lost I)", 
                      col = "number of \ngenerations \n(ngens)") + theme_classic() 
freqR_plot <- ggplot(tgr, aes(unlist(reps), unlist(avg_freq_R_tf), col = unlist(gens))) + 
  geom_point() + labs(x = "number of repetitions (bnum)", 
                      y = "freq R", # note error in code... didn't divide by 2 in original model
                      col = "number of \ngenerations \n(ngens)") + theme_classic() 

jpeg("plots/convg_check_extinct.jpeg",height=170,width=170,res=400,units='mm')
print(extinct_plot)
dev.off()

jpeg("plots/convg_check_timeK.jpeg",height=170,width=170,res=400,units='mm')
print(timeK_plot) # why is it splitting into two based on number of gens?
dev.off()

jpeg("plots/convg_check_lostI.jpeg",height=170,width=170,res=400,units='mm')
print(lostI_plot)
dev.off()

jpeg("plots/convg_check_freqR_tf.jpeg",height=170,width=170,res=400,units='mm')
print(freqR_plot)
dev.off()

