library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(patchwork)
library(emmeans)

setwd("/Volumes/Seagate Por")

surv_loc_df=read.table("PhD/Chapter_5_spatial_ID_x_E/Spatial_var_inbreeding/survival_loc.txt", sep = ",", header = TRUE)%>%
  select(-MumFROH, -BirthWt)%>%
  na.omit()

head(surv_loc_df)

surv_loc_df$Code=as.factor(surv_loc_df$Code)
surv_loc_df$MumCode=as.factor(surv_loc_df$MumCode)
surv_loc_df$BirthYear=as.factor(surv_loc_df$BirthYear)
surv_loc_df$Sex=as.factor(surv_loc_df$Sex)
surv_loc_df$MotherStatus=as.factor(surv_loc_df$MotherStatus)
surv_loc_df$Reg=as.factor(surv_loc_df$Reg)
## base model of juvenile survival
suv_model_simple=glmmTMB(juvenile_survival~ 1+ Sex + MotherStatus + mum_age+mum_age_sq+Day_seq+FROH+
                           (1|BirthYear)+(1|MumCode), 
                         family=binomial, 
                         data=surv_loc_df, 
                         na.action = na.omit)


surv_reg_fixed=update(suv_model_simple, ~ . + Reg) ##just region as fixed effect
summary(surv_reg_fixed)


## predictions of juvenile survival propability across the 6 regions 
# between sexes 
reg_pred_f=ggpredict(surv_reg_fixed, terms = c("Reg", "Sex[1]"))%>%
  mutate(x = fct_relevel(x, "SI", "IM", "LA", "NG","MG",  "SG"))
reg_pred_m=ggpredict(surv_reg_fixed, terms = c("Reg", "Sex[2]"))%>%
  mutate(x = fct_relevel(x, "SI", "IM", "LA", "NG","MG",  "SG"))
reg_pred=rbind(reg_pred_f,reg_pred_m)%>%
  arrange(x)%>%
  na.omit()

surv_reg=reg_pred%>%
  ggplot(aes(x=x, y=predicted, color=x, ymin=conf.low, ymax=conf.high, group=group))+
  geom_pointrange(linewidth=1, position = position_dodge(width=0.5))+
  theme_bw()+
  scale_color_manual(values = c("#f0a30a" ,"#a20025","#00aba9","chocolate1","#60a917", "#647687"))+
  labs(x="Spatial region", y="Predicted juvenile survival probability")+
  theme(text = element_text(size = 18),legend.position = "none")
surv_reg


##################################
## INTERACTION OF F WITH REGION ##
##################################

surv_froh_inter=update(suv_model_simple, ~ . -FROH + (Reg*FROH)) #interaction between region and froh
summary(surv_froh_inter)

#plot predictions
inter=plot(predict_response(surv_froh_inter, terms = c("FROH[all]","Reg")),show.title=FALSE, line.size=1, colors="metro")+
  labs(x = expression(F["ROH"]), y = "Juvenile survival probability", colour = "Spatial \nregion")+
  theme(text = element_text(size = 15)) +
  xlim(0,0.2)
inter



# only testing slope difference between 0-0.2 FROH because below this is where there are very few ids 
test_predictions(surv_froh_inter,c("FROH[0,0.05,0.1,0.15,0.2]","Reg"), p_adjust = "bonferroni")

# 
# Reg   | Contrast |       95% CI |      p
# ----------------------------------------
#   IM-LA |    -0.75 | -2.65,  1.16 | > .999
# IM-MG |    -1.93 | -3.63, -0.24 | 0.378 
# IM-NG |    -0.54 | -2.02,  0.95 | > .999
# IM-SG |    -2.76 | -4.34, -1.18 | 0.009  #### <<<<<<
# IM-SI |    -0.62 | -2.24,  1.01 | > .999
# LA-MG |    -1.19 | -3.27,  0.90 | > .999
# LA-NG |     0.21 | -1.72,  2.14 | > .999
# LA-SG |    -2.01 | -4.01, -0.01 | 0.729 
# LA-SI |     0.13 | -1.92,  2.18 | > .999
# MG-NG |     1.40 | -0.32,  3.12 | > .999
# MG-SG |    -0.82 | -2.61,  0.97 | > .999
# MG-SI |     1.32 | -0.54,  3.18 | > .999
# NG-SG |    -2.22 | -3.83, -0.61 | 0.103 
# NG-SI |    -0.08 | -1.74,  1.58 | > .999
# SG-SI |     2.14 |  0.40,  3.88 | 0.243 





