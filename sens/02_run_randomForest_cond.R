library(dplyr)
library(randomForest)
library(ggplot2)
library(cowplot)

# loading the data -- note weird directory because of HPC use
sum_0302a <- readRDS("dat/sim_dat_0302.Rdata")
sum_0302b <- readRDS("dat/sim_dat_2_0302.Rdata")
random_parms <- readRDS("dat/mat_var_0211.Rdata") # figure out what these values were converging around....? new values have super little er
# ^ read correct matrix!

# convert to df
sim_dat <- rbind(sum_0302a, sum_0302b) # just one???
# filter out P>50 whole time
sim_dat_drop <- sim_dat %>% filter(pop_drop50==1) %>% 
  group_by(parm_number) %>% mutate(rep_num = n())
# calc p(ex), p(er) & time, p(il)
sim_dat_drop_ext <- sim_dat_drop %>%
  filter(extinct == 1) %>% 
  group_by(parm_number) %>% 
  summarise(`P(extinct)` = n()/rep_num[1])
sim_dat_drop_er <- sim_dat_drop %>%
  filter(extinct == 0, at_K95 == 1, r_allele_peak45 == 1) %>% 
  group_by(parm_number) %>% 
  # already filtered pop_drop50
  summarise(`P(ER)` = n()/rep_num[1], ER_K95 = mean(firstK95), ER_IF = mean(final_inf_prev))
sim_dat_drop_il <- sim_dat_drop %>%
  filter(extinct == 0, at_K95 == 1, r_allele_peak45 == 0, final_inf_prev < 0.15) %>% 
  # already filtered pop_drop50, using generous final_inf_prev to capture declining disease, e.g., for SIR declines
  group_by(parm_number) %>% 
  summarise(`P(IL)` = n()/rep_num[1], IL_IF = mean(final_inf_prev))
sim_dat_drop_p50 <- sim_dat_drop %>% 
  group_by(parm_number) %>% 
  summarise(`P(drop)` = rep_num[1]/1000)
# merging things
tmp <- merge(sim_dat_drop_ext, sim_dat_drop_er, by = "parm_number", all = T)
tmp <- merge(tmp, sim_dat_drop_il, by = "parm_number", all = T)
tmp <- merge(tmp, sim_dat_drop_p50, by = "parm_number", all = T)
tmp$tot_probs <- apply(tmp[, c(2, 3, 6)], 1, sum, na.rm = T)
GSA_cond_dat <- tmp 

# remove NAs for P() --> those are 0s
GSA_cond_dat$`P(extinct)`[is.na(GSA_cond_dat$`P(extinct)`)] <- 0
GSA_cond_dat$`P(ER)`[is.na(GSA_cond_dat$`P(ER)`)] <- 0
GSA_cond_dat$`P(IL)`[is.na(GSA_cond_dat$`P(IL)`)] <- 0
GSA_cond_dat$unkn <- 1-GSA_cond_dat$tot_probs

# load in the parameter files
RF_dat <- merge(as.data.frame(random_parms), GSA_cond_dat, by.x = "V36", by.y = "parm_number")
# rename so it's easy to navigate and plot
colnames(RF_dat) <- c("parameter number", "compartments", "event order", 
                      "\u03B2_F,RR", "\u03B2_F,WR", "\u03B2_F,WW", "\u03C3_\u03B2,F", 
                      "\u03B2_D,RR", "\u03B2_D,WR", "\u03B2_D,WW", "\u03C3_\u03B2,D", 
                      "\u03BC_S,RR", "\u03BC_S,WR", "\u03BC_S,WW", "\u03BC_I,RR", "\u03BC_I,WR", "\u03BC_I,WW", "\u03C3_\u03BC,S", "\u03C3_\u03BC,I", 
                      "\u03B3_RR", "\u03B3_WR", "\u03B3_WW", "\u03C3_\u03B3", 
                      "\u03BB_RR", "\u03BB_WR", "\u03BB_WW", "\u03BB_sd", "\u03BD",
                      "K", "\u03C3_K",
                      "t_q", "d", 
                      "Q_0", "I_0", "t_f", "transmission case", colnames(RF_dat[37:45]))
# note the first two variables are categorical
RF_dat$"compartments" <- as.factor(RF_dat$"compartments"); RF_dat$"event order" <- as.factor(RF_dat$"event order"); RF_dat$"transmission case" <- as.factor(RF_dat$"transmission case")

# trying with ratios instead of raw values, e.g., ratio of RR:WW
RF_ratio <- data.frame(compartments = RF_dat$"compartments", 
                       "event order" = RF_dat$"event order", 
                       "allele benefit, enviro. trans. blocking" = RF_dat$"\u03B2_F,RR"/RF_dat$"\u03B2_F,WW", 
                       "wild-type enviro. tranmission" = RF_dat$"\u03B2_F,WW", 
                       "enviro. tranmission variation" = RF_dat$"\u03C3_\u03B2,F",
                       "transmission type" = RF_dat$"transmission case", # want to be able to seperate cases with and without enviro
                       "allele benefit, dens. trans. blocking" = RF_dat$"\u03B2_D,RR"/RF_dat$"\u03B2_D,WW",
                       "wild-type dens. tranmission" = RF_dat$"\u03B2_D,WW",
                       "dens tranmission variation" = RF_dat$"\u03C3_\u03B2,D",
                       "wild-type natural mortality" = RF_dat$"\u03BC_S,WW",
                       "allele benefit, mortality-blocking" = RF_dat$"\u03BC_I,RR"/RF_dat$"\u03BC_I,WW",
                       "wild-type disease mortality" = RF_dat$"\u03BC_I,WW",
                       "natural mortality variation" = RF_dat$"\u03C3_\u03BC,S", 
                       "disease mortality variation" = RF_dat$"\u03C3_\u03BC,I", 
                       "allele benefit, clearance aug." = RF_dat$"\u03B3_RR"/RF_dat$"\u03B3_WW", 
                       "wild-type recovery" = RF_dat$"\u03B3_WW", 
                       "recovery variation" = RF_dat$"\u03C3_\u03B3", 
                       "allele cost, fecundity" = RF_dat$"\u03BB_RR"/RF_dat$"\u03BB_WW",
                       "wild-type fecundity" = RF_dat$"\u03BB_WW"
) # then everything after is the same
RF_ratio <- cbind(RF_ratio, RF_dat[, 27:35]) # , RF_dat[,c(43:44)]) # RF_dat add is for and r/inf info
colnames(RF_ratio) <- c("compartments", "event order",
                        "allele benefit: env. trans. block.", "wild type env. tranmission", "env. tranmission variation", 
                        "transmission type", 
                        "allele benefit: dens. trans. block.", "wild type dens. tranmission", "dens. tranmission variation",
                        "wild type natural mortality", "allele benefit: mortality block.", "wild type disease mortality", "natural mortality variation", "disease mortality variation",
                        "allele benefit: clearance aug.", "wild type recovery", "recovery variation",
                        "allele cost: fecundity", "wild type fecundity", "fecundity variation", 
                        "mutation rate", 
                        "carrying capacity", "carrying capactity variation",
                        "initalization length",
                        "disease gens. per host gen.", 
                        "allele frequency after initalization", "initial disease prevelence", "total simulation time")
RF_disc <- data.frame(name = colnames(RF_ratio), 
                      typ = c("disease ecology", "force of disease", 
                              "adaptation", "force of disease", "other", 
                              "disease ecology",
                              "adaptation", "force of disease", "other", 
                              "other", "adaptation", "force of disease", "other", "other",
                              "adaptation", "force of disease", "other", 
                              "adaptation", "other", "other", 
                              "adaptation", 
                              "other", "other", "other", 
                              "force of disease", "other", "other", "other"))
saveRDS(RF_disc, file = "dat/RF_names.Rdata")

set.seed(108) # for consistency in RF calcs

# extinction
EX_cond <- randomForest(x=RF_ratio, y=RF_dat$`P(extinct)`, ntree=1600, importance = T, localImp = T)
plot(EX_cond) # this is checking the convergence?
varImpPlot(EX_cond)
EX_cond # checing var explained etc
saveRDS(EX_cond, file = "dat/EX_RF_cond.Rdata")

# ER
ER_cond <- randomForest(x=RF_ratio, y=RF_dat$`P(ER)`, ntree=1600, importance = T, localImp = T)
plot(ER_cond) # this is checking the convergence?
varImpPlot(ER_cond)
ER_cond # checing var explained etc
saveRDS(ER_cond, file = "dat/ER_RF_cond.Rdata")

# IL
IL_cond <- randomForest(x=RF_ratio, y=RF_dat$`P(IL)`, ntree=1600, importance = T, localImp = T)
plot(IL_cond) # this is checking the convergence?
varImpPlot(IL_cond)
IL_cond # checing var explained etc
saveRDS(IL_cond, file = "dat/IL_RF_cond.Rdata")

# UNK_cond <- randomForest(x=RF_ratio, y=RF_dat$unkn, ntree=1600, importance = T, localImp = T)
# plot(UNK_cond) # this is checking the convergence?
# varImpPlot(UNK_cond)
# UNK_cond # checing var explained etc
# saveRDS(UNK_cond, file = "dat/UNK_RF_cond.Rdata")

# and now formal plots -- take the top 8
ex_dat <- as.data.frame(importance(EX_cond)); ex_dat$names <- RF_disc$name; ex_dat$typ <- RF_disc$typ
ex_dat$case <- rep(paste0("P(Ext), ", round(EX_cond$rsq[1000]*100, 2), "% var exp"))
ex_dat <- ex_dat %>% arrange(desc(`%IncMSE`)); ex_dat <- ex_dat[1:10, ]
EX_RF_plt <- ggplot(data = ex_dat, aes(`%IncMSE`, reorder(names, `%IncMSE`))) + 
  geom_linerange(aes(xmin = 0, xmax = `%IncMSE`)) + 
  geom_point(aes(col = typ), size = 3) + # alt: aes(size = log(IncNodePurity))
  labs(x = "", y = NULL, col = "parameter\ntype:") + 
  scale_color_manual(values = c("#ac1457", "#DB6341", "#f1c4a2"), 
                     breaks = c("adaptation", "force of disease", "disease ecology"), 
                     labels = c("adaptation", "force\nof disease", "disease\necology")) + 
  theme_bw() + facet_wrap(~case) +  
  theme(text = element_text(size = 10), 
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8), 
        legend.position = "bottom") 

il_dat <- as.data.frame(importance(IL_cond)); il_dat$names <- RF_disc$name; il_dat$typ <- RF_disc$typ
il_dat$case <- rep(paste0("P(IL), ", round(IL_cond$rsq[1000]*100, 2), "% var exp"))
il_dat <- il_dat %>% arrange(desc(`%IncMSE`)); il_dat <- il_dat[1:10, ]
IL_RF_plt <- ggplot(data = il_dat, aes(`%IncMSE`, reorder(names, `%IncMSE`))) + 
  geom_linerange(aes(xmin = 0, xmax = `%IncMSE`)) + 
  geom_point(aes(col = typ), size = 3) + # alt: aes(size = log(IncNodePurity))
  labs(x = "", y = NULL, col = "parameter type:") + 
  scale_color_manual(values = c("#ac1457", "#DB6341", "#f1c4a2"), breaks = c("adaptation", "force of disease", "disease ecology")) + 
  theme_bw() + facet_wrap(~case) +  
  theme(text = element_text(size = 10)) 

er_dat <- as.data.frame(importance(ER_cond)); er_dat$names <- RF_disc$name; er_dat$typ <- RF_disc$typ
er_dat$case <- rep(paste0("P(ER), ", round(ER_cond$rsq[1000]*100, 2), "% var exp"))
er_dat <- er_dat %>% arrange(desc(`%IncMSE`)); er_dat <- er_dat[1:10, ]
ER_RF_plt <- ggplot(data = er_dat, aes(`%IncMSE`, reorder(names, `%IncMSE`))) + 
  geom_linerange(aes(xmin = 0, xmax = `%IncMSE`)) + 
  geom_point(aes(col = typ), size = 3) + # alt: aes(size = log(IncNodePurity))
  labs(x = "Importance", y = NULL, col = "parameter type:") + 
  scale_color_manual(values = c("#ac1457", "#DB6341", "#f1c4a2"), breaks = c("adaptation", "force of disease", "disease ecology")) + 
  theme_bw() + facet_wrap(~case) +  
  theme(text = element_text(size = 10)) 

legend <- get_plot_component(EX_RF_plt, 'guide-box-bottom', return_all = TRUE)
RF_plt <- plot_grid(EX_RF_plt + theme(legend.position = "none"), 
                     IL_RF_plt + theme(legend.position = "none"),
                     ER_RF_plt + theme(legend.position = "none"),
                     legend,
                     ncol = 1, rel_heights = c(1, 1, 1, 0.15))
# save figure
png("figs/figure_plot/RF_cond.png",height=200,width=85,res=400,units='mm')
print(RF_plt)
dev.off()

