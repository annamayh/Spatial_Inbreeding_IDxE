library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(patchwork)
library(emmeans)


setwd("/Volumes/Seagate_HD/")

bw_df=read.table("Deer_spatial_variation_ID/birthweight_loc_2024_new_calves.txt", sep = ",", header = TRUE)%>%
  drop_na(N, E)%>%
  mutate(N_scaled=N/100,
         E_scaled=E/100,
         FROH_trans=sqrt(FROH))

bw_df$Code=as.factor(bw_df$Code)
bw_df$MumCode=as.factor(bw_df$MumCode)
bw_df$BirthYear=as.factor(bw_df$BirthYear)
bw_df$Sex=as.factor(bw_df$Sex)
bw_df$Reg=as.factor(bw_df$Reg)

bw_df$N_scaled=as.numeric(bw_df$N_scaled)
bw_df$E_scaled=as.numeric(bw_df$E_scaled)




bw_model_simple=glmmTMB(CaptureWt~ Sex + AgeHrs+ MotherStatus+mum_age+mum_age_sq+Day_seq+
                          (1|BirthYear)+ (1|MumCode), 
                        family=gaussian(), 
                        data=bw_df, 
                        na.action = na.omit,
)


bw_froh_inter=update(bw_model_simple, ~ . + (N_scaled*FROH_trans)+(E_scaled*FROH_trans)) #interaction between region and froh
summary(bw_froh_inter)

quantile(bw_df$N, na.rm = T)
sd(bw_df$N)
mean(bw_df$N)

inter_bw_N=predict_response(bw_froh_inter, terms = c("FROH_trans[all]","N_scaled[80.44, 80.23]"), type = 'fixed')

inter_bw_N_plot=inter_bw_N%>%
  as.data.frame()%>%
  mutate(x_rescaled = x^2, 
         group=as.numeric(as.character(group)),
         northing_rescaled = as.factor(group*100))%>%
  ggplot(aes(x=x_rescaled, y=predicted, ymin=conf.low, ymax=conf.high,
             group = group, colour = northing_rescaled, fill=northing_rescaled))+
  theme_bw()+
  geom_line(line_width=1, linetype='dashed', show.legend = F)+
  geom_ribbon(alpha=0.4, colour = NA, show.legend = F)+
  labs(x = expression(F["ROH"]), y = "Predicted birth weight (kg)", colour = "Northing")+
  theme(text = element_text(size = 15)) +
  scale_colour_discrete(labels = c("Mean - SD", "Mean + SD"))


inter_bw_N_plot




save(inter_bw_N_plot, file = "Deer_spatial_variation_ID/plots/bw_IDxN.RData")

