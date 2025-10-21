library(dplyr)
library(randomForest)
library(ggplot2)
library(cowplot)

# UNK_cond <- randomForest(x=RF_ratio, y=RF_dat$unkn, ntree=1600, importance = T, localImp = T)
# plot(UNK_cond) # this is checking the convergence?
# varImpPlot(UNK_cond)
# UNK_cond # checing var explained etc
# saveRDS(UNK_cond, file = "dat/UNK_RF_cond.Rdata")

# load in randomforests
EX <- readRDS("dat/RF_ext_ratio.Rdata") # etc.
IL <- forest_DR_ratio
ER <- forest_ER_ratio

RF_class <- data.frame(name_vals = colnames(RF_ratio), 
                       typ = c("disease ecology", "disease ecology", "force of disease", "disease ecology",
                               "force of disease", "background", "force of disease", "force of disease", 
                               "adaptation", "background", "adaptation", "adaptation", "adaptation", 
                               "background", "adaptation", "background", "background", "disease ecology", "background"))

# and now formal plots -- take the top 8
ex_dat <- as.data.frame(importance(EX)); ex_dat$name_vals <- RF_class$name_vals; ex_dat$typ <- RF_class$typ
ex_dat$case <- rep(paste0("P(Ext), ", round(EX$rsq[1000]*100, 2), "% var exp"))
ex_dat <- ex_dat %>% arrange(desc(`%IncMSE`)); ex_dat <- ex_dat[1:10, ]
EX_RF_plt <- ggplot(data = ex_dat, aes(`%IncMSE`, reorder(name_vals, `%IncMSE`))) + 
  geom_linerange(aes(xmin = 0, xmax = `%IncMSE`)) + 
  geom_point(aes(col = typ), size = 3) + # alt: aes(size = log(IncNodePurity))
  labs(x = "", y = NULL, col = "parameter\ntype:") + 
  scale_color_manual(values = c("#ac1457", "#DB6341", "#f1c4a2", "black"), 
                     breaks = c("adaptation", "force of disease", "disease ecology", "background"), 
                     labels = c("adaptation", "force\nof disease", "disease\necology", "background")) + 
  theme_bw() + facet_wrap(~case) +  
  theme(text = element_text(size = 12), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12), 
        legend.position = "bottom") 

il_dat <- as.data.frame(importance(IL)); il_dat$name_vals <- RF_class$name_vals; il_dat$typ <- RF_class$typ
il_dat$case <- rep(paste0("P(IL), ", round(IL$rsq[1000]*100, 2), "% var exp"))
il_dat <- il_dat %>% arrange(desc(`%IncMSE`)); il_dat <- il_dat[1:10, ]
IL_RF_plt <- ggplot(data = il_dat, aes(`%IncMSE`, reorder(name_vals, `%IncMSE`))) + 
  geom_linerange(aes(xmin = 0, xmax = `%IncMSE`)) + 
  geom_point(aes(col = typ), size = 3) + # alt: aes(size = log(IncNodePurity))
  labs(x = "", y = NULL, col = "parameter type:") + 
  scale_color_manual(values = c("#ac1457", "#DB6341", "#f1c4a2", "black"), 
                     breaks = c("adaptation", "force of disease", "disease ecology", "background"), 
                     labels = c("adaptation", "force\nof disease", "disease\necology", "background")) + 
  theme_bw() + facet_wrap(~case) +  
  theme(text = element_text(size = 10)) 

er_dat <- as.data.frame(importance(ER)); er_dat$name_vals <- RF_class$name_vals; er_dat$typ <- RF_class$typ
er_dat$case <- rep(paste0("P(ER), ", round(ER$rsq[1000]*100, 2), "% var exp"))
er_dat <- er_dat %>% arrange(desc(`%IncMSE`)); er_dat <- er_dat[1:10, ]
ER_RF_plt <- ggplot(data = er_dat, aes(`%IncMSE`, reorder(name_vals, `%IncMSE`))) + 
  geom_linerange(aes(xmin = 0, xmax = `%IncMSE`)) + 
  geom_point(aes(col = typ), size = 3) + # alt: aes(size = log(IncNodePurity))
  labs(x = "Importance", y = NULL, col = "parameter type:") + 
  scale_color_manual(values = c("#ac1457", "#DB6341", "#f1c4a2", "black"), 
                     breaks = c("adaptation", "force of disease", "disease ecology", "background"), 
                     labels = c("adaptation", "force\nof disease", "disease\necology", "background")) + 
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

