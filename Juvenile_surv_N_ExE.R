library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(patchwork)
library(emmeans)
library(see)

setwd("/Volumes/Seagate_HD/")

surv_loc_df=read.table("Deer_spatial_variation_ID/survival_loc_2024_new_calves.txt", sep = ",", header = TRUE)%>%
  drop_na(N, E)%>%
  mutate(N_scaled=N/100,
         E_scaled=E/100,
         FROH_trans=sqrt(FROH))%>%
  select(-MumFROH)
head(surv_loc_df)

surv_loc_df%>%
  group_by(Reg)%>%
  summarise(mean=mean(juvenile_survival))


surv_loc_df$Code=as.factor(surv_loc_df$Code)
surv_loc_df$MumCode=as.factor(surv_loc_df$MumCode)
surv_loc_df$BirthYear=as.factor(surv_loc_df$BirthYear)
surv_loc_df$Sex=as.factor(surv_loc_df$Sex)
surv_loc_df$MotherStatus=as.factor(surv_loc_df$MotherStatus)
surv_loc_df$Reg=as.factor(surv_loc_df$Reg)

surv_loc_df$N_scaled=as.numeric(surv_loc_df$N_scaled)
surv_loc_df$E_scaled=as.numeric(surv_loc_df$E_scaled)


## base model of juvenile survival
suv_model_simple=glmmTMB(juvenile_survival~ 1+ Sex + MotherStatus + mum_age+mum_age_sq+Day_seq+
                           (1|BirthYear)+(1|MumCode), 
                         family=binomial(), 
                         data=surv_loc_df, 
                         na.action = na.omit)

summary(suv_model_simple)

surv_froh_N_E=update(suv_model_simple, ~ . + FROH_trans+ N_scaled + E_scaled) #interaction between region and froh
summary(surv_froh_N_E)

surv_froh_inter=update(suv_model_simple, ~ . + N_scaled + E_scaled + FROH_trans+(FROH_trans:N_scaled)+(FROH_trans:E_scaled)) #interaction between region and froh
summary(surv_froh_inter)

plot(ggpredict(surv_froh_inter, terms = c("N_scaled","FROH_trans")))# bascially shows that the 

inter_surv_N=ggemmeans(surv_froh_inter, terms = c("FROH_trans[all]","N_scaled[80.44, 80.23]"))
                  
inter_surv_N_plot=inter_surv_N%>%
  as.data.frame()%>%
  mutate(x_rescaled = x^2, 
         group=as.numeric(as.character(group)),
         northing_rescaled = as.factor(group*100))%>%
  ggplot(aes(x=x_rescaled, y=predicted, ymin=conf.low, ymax=conf.high,
                    group = northing_rescaled, colour = northing_rescaled, fill=northing_rescaled))+
  theme_bw()+
  geom_line(line_width=1)+
  geom_ribbon(alpha=0.4, colour = NA, show.legend = F)+
  labs(x = expression(F["ROH"]), y = "Predicted juvenile survival probability", colour = "Northing")+
  theme(text = element_text(size = 15), 
        plot.title = element_text(hjust=0.5)) +
  scale_colour_discrete(labels = c("Mean - SD", "Mean + SD"))
                  
inter_surv_N_plot


# ggsave(inter_surv_N_plot, 
#        file = "Deer_spatial_variation_ID/REVISED/Figs/Fig.5.png",
#        width = 6,
#        height = 5)

both=inter_bw_N_plot+inter_surv_N_plot+plot_annotation(tag_levels = 'A')+
  plot_layout(guides = "collect")
both

## get per region
surv_loc_df%>%
  group_by(Reg)%>%
  summarise(mean=mean(N))




ggsave(both,
       file = "Deer_spatial_variation_ID/REVISED/IDxN_SUPP.png",
       width = 10,
       height = 5)
