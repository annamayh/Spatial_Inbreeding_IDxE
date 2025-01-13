library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(patchwork)
library(emmeans)

setwd("/Volumes/Seagate_HD/")

surv_loc_df=read.table("Deer_spatial_variation_ID/survival_loc_2024_new_calves.txt", sep = ",", header = TRUE)%>%
  mutate(FROH_trans=sqrt(FROH))

table(surv_loc_df$juvenile_survival)

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

## base model of juvenile survival
suv_model_simple=glmmTMB(juvenile_survival~ 1+ Sex + MotherStatus + mum_age+mum_age_sq+Day_seq+FROH_trans+
                           (1|BirthYear)+(1|MumCode), 
                         family=binomial(), 
                         data=surv_loc_df, 
                         na.action = na.omit)


surv_reg_fixed=update(suv_model_simple, ~ . + Reg) ##just region as fixed effect
summary(surv_reg_fixed)

##################################
## INTERACTION OF F WITH REGION ##
##################################

surv_froh_inter=update(suv_model_simple, ~ . -FROH_trans + (Reg*FROH_trans)) #interaction between region and froh
summary(surv_froh_inter)

  

## emmtrends 
emmeans::emtrends(surv_froh_inter, pairwise ~ Reg, var = "FROH_trans")


pred_surv=predict_response(surv_froh_inter, terms = c("FROH_trans[all]","Reg"))%>%
  as.data.frame()%>%
  mutate(FROH=x^2)## tranforming back into FROH

cols = c("SI"="#f0a30a" ,"IM" = "#a20025", "LA" = "#1ba1e2", "NG" = "#fa6800", "MG" = "#008a00", "SG" = "#647687")

inter_surv=ggplot(pred_surv, aes(x=FROH, y=predicted, ymin=conf.low, ymax=conf.high,
                                                   group = group, colour = group, fill=group))+
  geom_line(linewidth=1)+
  geom_ribbon(alpha=0.15, colour = NA, show.legend = F)+
  labs(x = expression(F["ROH"]), y = "Predicted juvenile survival probability", colour = "Spatial \nregion")+
  theme_bw()+
  theme(text = element_text(size = 15)) +
  scale_color_manual(values = cols, breaks=c("SI","IM","LA", "NG", "MG", "SG"))+
  scale_fill_manual(values = cols)

inter_surv





#test_predictions(surv_froh_inter,c("FROH[0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325]","Reg"), p_adjust = "bonferroni")



# ggsave(inter,
#        file = "Deer_spatial_variation_ID/plots/IDxE_surv.png",
#        width = 7,
#        height = 6, 
#        bg = "white"
#       )
# 
# 
# ggsave(inter_extended,
#        file = "Deer_spatial_variation_ID/plots/IDxE_surv_SUPP.png",
#        width = 7,
#        height = 6, 
#        bg = "white"
# )

save(inter_surv, file = "Deer_spatial_variation_ID/plots/surv_IDxE_plots.RData")
# updated 25.11.23
