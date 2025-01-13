library(tidyverse)
library(patchwork)

setwd("/Volumes/Seagate_HD/")

## FIG 3 ##
load("Deer_spatial_variation_ID/plots/surv_plots.RData")
load("Deer_spatial_variation_ID/plots/bw_plots.RData")


all=(bw_reg+inla_bw_plot)/(surv_reg+inla_surv_plot)+plot_annotation(tag_levels = 'A')
all

ggsave(all,
       file = "Deer_spatial_variation_ID/REVISED/Figs/Fig.3.jpeg",
       width = 11,
       height = 12,
       dpi=300)



## FIG 4 ##
load("Deer_spatial_variation_ID/plots/bw_IDxE.RData")
load("Deer_spatial_variation_ID/plots/surv_IDxE_plots.RData")

surv_loc_df=read.table("Deer_spatial_variation_ID/survival_loc_2024_new_calves.txt", sep = ",", header = TRUE)%>%
  filter(!E>1385)%>%
  filter(!N<7997.5)#remove 38 records outside the syudy areas


## FROH per region violin plots 

cols <- c("SI"="#f0a30a" ,"IM"="#a20025","LA"="#00aba9","NG"="chocolate1","MG"="#60a917", "SG"="#647687")

violin=ggplot(surv_loc_df, aes(x=Reg, y=FROH))+
  scale_x_discrete(limits=c("SG", "MG", "NG", "LA", "IM", "SI"))+
  geom_jitter(aes(colour=Reg), position=position_jitter(0.2), alpha=0.5)+
  scale_color_manual(values=cols)+
  geom_violin(fill = "transparent")+
  coord_flip()+
  theme_bw()+
  geom_hline(yintercept = 0.1, linetype=2)+
  geom_hline(yintercept = 0.2, linetype=2)+
  geom_hline(yintercept = 0.3, linetype=2)+
  theme(text = element_text(size = 14), 
        axis.title.y = element_blank(),
        #axis.text.y = element_blank(),
        legend.position = "none")+
  labs(y=expression(F["ROH"])) +
  stat_summary(fun = "mean",
               geom = "point",
               color = "black")

violin

IDxE=(inter_bw+inter_surv)/violin+
  plot_annotation(tag_levels = 'A')+plot_layout(guides = "collect") +
  plot_layout(heights = c(7,4))
IDxE

ggsave(IDxE,
       file = "Deer_spatial_variation_ID/REVISED/Figs/Fig.4.jpeg",
       width = 9,
       height = 8,
       dpi=300)


