library(ggplot2)
library(dplyr)
library(cowplot)
load("figs/sim_fig_dat/forest_ER_ratio")
load("figs/sim_fig_dat/forest_extinct_ratio")
load("figs/sim_fig_dat/forest_K95_ratio")
# also need to run the RF_ratio in 02_run_randomForest.R

# and now formal plots -- take the top 8
ex_dat <- as.data.frame(importance(forest_extinct_ratio)); ex_dat$names <- colnames(RF_ratio)
ex_dat <- ex_dat %>% arrange(desc(`%IncMSE`)); ex_dat <- ex_dat[1:8, ]
RF_ex <- ggplot(data = ex_dat, aes(`%IncMSE`, reorder(names, `%IncMSE`))) + 
  geom_linerange(aes(xmin = 0, xmax = `%IncMSE`)) + 
  geom_point(col = "#ac1457", size = 3) + # alt: aes(size = log(IncNodePurity))
  labs(x = "", y = NULL, title = "P(Extinct), 76.04%") + 
  theme_bw()  + 
  theme(text = element_text(size = 10), 
        plot.title = element_text(hjust = 0.5, size = 10)) 

er_dat <- as.data.frame(importance(forest_ER_ratio)); er_dat$names <- colnames(RF_ratio)
er_dat <- er_dat %>% arrange(desc(`%IncMSE`)); er_dat <- er_dat[1:8, ]
RF_er <- ggplot(data = er_dat, aes(`%IncMSE`, reorder(names, `%IncMSE`))) + 
  geom_linerange(aes(xmin = 0, xmax = `%IncMSE`)) + 
  geom_point(col = "#DB6341", size = 3) + 
  labs(x = "", y = NULL, title = "P(Evolutionary Rescue), 31.84%") + 
  theme_bw()  + 
  theme(text = element_text(size = 10), 
        plot.title = element_text(hjust = 0.5, size = 10)) 

K95_dat <- as.data.frame(importance(forest_K95_ratio)); K95_dat$names <- colnames(RF_ratio)
K95_dat <- K95_dat %>% arrange(desc(`%IncMSE`)); K95_dat <- K95_dat[1:8, ]
RF_K95 <- ggplot(data = K95_dat, aes(`%IncMSE`, reorder(names, `%IncMSE`))) + 
  geom_linerange(aes(xmin = 0, xmax = `%IncMSE`)) + 
  geom_point(col = "#f1c4a2", size = 3) + 
  labs(x = "Importance", y = NULL, title = "Time to rescue, 40.82%") + 
  theme_bw()  + 
  theme(text = element_text(size = 10), 
        plot.title = element_text(hjust = 0.5, size = 10)) 

RF_plot <- plot_grid(RF_ex, RF_er, RF_K95, ncol = 1)
# save figure
png("figs/figure_plot/random_forest.png",height=155,width=85,res=400,units='mm')
print(RF_plot)
dev.off()

# note option to see one at a time effects  (no interactions)
partialPlot(forest_extinct_ratio, RF_ratio, "λ_ratio") 
