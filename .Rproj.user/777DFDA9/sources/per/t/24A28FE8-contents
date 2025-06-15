library(ggplot2)
library(dplyr)
library(cowplot)
library(randomForest)
RF_ER <- readRDS("dat/RF_ER_ratio.Rdata")
RF_K95 <- readRDS("dat/RF_K95_ratio.Rdata")
RF_IL <- readRDS("dat/RF_IL_ratio.Rdata")
RF_ext <- readRDS("dat/RF_ext_ratio.Rdata")
RF_names <- readRDS("dat/RF_names.Rdata")
# load("dat/forest_extinct_ratio")
# also need to run the RF_ratio in 02_run_randomForest.R to get col names

# and now formal plots -- take the top 10
ex_dat <- as.data.frame(importance(RF_ext)); ex_dat$names <- RF_names
ex_dat$case <- rep(paste0("P(Ext.), ", round(RF_ext$rsq[1000]*100, 2), "% var exp"))
ex_dat <- ex_dat %>% arrange(desc(`%IncMSE`)); ex_dat <- ex_dat[1:8, ]
EX_RF_plt <- ggplot(data = ex_dat, aes(`%IncMSE`, reorder(names, `%IncMSE`))) + 
  geom_linerange(aes(xmin = 0, xmax = `%IncMSE`)) + 
  geom_point(col = "black", size = 3) + # alt: aes(size = log(IncNodePurity))
  labs(x = "", y = NULL) + 
  theme_bw() + facet_wrap(~case) +  
  theme(text = element_text(size = 8), 
        plot.title = element_text(hjust = 0.5, size = 8)) 

il_dat <- as.data.frame(importance(RF_IL)); il_dat$names <- RF_names
il_dat$case <- rep(paste0("P(IL), ", round(RF_IL$rsq[1000]*100, 2), "% var exp"))
il_dat <- il_dat %>% arrange(desc(`%IncMSE`)); il_dat <- il_dat[1:8, ]
IL_RF_plt <- ggplot(data = il_dat, aes(`%IncMSE`, reorder(names, `%IncMSE`))) + 
  geom_linerange(aes(xmin = 0, xmax = `%IncMSE`)) + 
  geom_point(col = "#ac1457", size = 3) + # alt: aes(size = log(IncNodePurity))
  labs(x = "", y = NULL) + 
  theme_bw() + facet_wrap(~case) + 
  theme(text = element_text(size = 8), 
        plot.title = element_text(hjust = 0.5, size = 8)) 

er_dat <- as.data.frame(importance(RF_ER)); er_dat$names <- RF_names
er_dat$case <- rep(paste0("P(ER), ", round(RF_ER$rsq[1000]*100, 2), "% var exp"))
er_dat <- er_dat %>% arrange(desc(`%IncMSE`)); er_dat <- er_dat[1:8, ]
ER_RF_plt <- ggplot(data = er_dat, aes(`%IncMSE`, reorder(names, `%IncMSE`))) + 
  geom_linerange(aes(xmin = 0, xmax = `%IncMSE`)) + 
  geom_point(col = "#DB6341", size = 3) + 
  labs(x = "", y = NULL) + 
  theme_bw() + facet_wrap(~case) +  
  theme(text = element_text(size = 8), 
        plot.title = element_text(hjust = 0.5, size = 8)) 

K95_dat <- as.data.frame(importance(RF_K95)); K95_dat$names <- RF_names
K95_dat$case <- rep(paste0("Time to ER, ", round(RF_K95$rsq[1000]*100, 2), "% var exp"))
K95_dat <- K95_dat %>% arrange(desc(`%IncMSE`)); K95_dat <- K95_dat[1:8, ]
K95_RF_plt <- ggplot(data = K95_dat, aes(`%IncMSE`, reorder(names, `%IncMSE`))) + 
  geom_linerange(aes(xmin = 0, xmax = `%IncMSE`)) + 
  geom_point(col = "#f1c4a2", size = 3) + 
  labs(x = "Importance", y = NULL) + 
  theme_bw() + facet_wrap(~case) + 
  theme(text = element_text(size = 8), 
        plot.title = element_text(hjust = 0.5, size = 8)) 

RF_plot <- plot_grid(EX_RF_plt, IL_RF_plt, ER_RF_plt, K95_RF_plt, ncol = 1)
# save figure
png("figs/figure_plot/random_forest.png",height=220,width=85,res=400,units='mm')
print(RF_plot)
dev.off()

# note option to see one at a time effects  (no interactions)
partialPlot(RF_ER, RF_ratio, "λ_ratio") 
