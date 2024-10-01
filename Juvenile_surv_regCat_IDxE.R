library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(patchwork)
library(emmeans)

setwd("/Volumes/Seagate_HD/")

surv_loc_df=read.table("Deer_spatial_variation_ID/survival_loc_2024.txt", sep = ",", header = TRUE)
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
suv_model_simple=glmmTMB(juvenile_survival~ 1+ Sex + MotherStatus + mum_age+mum_age_sq+Day_seq+FROH+
                           (1|BirthYear)+(1|MumCode), 
                         family=binomial(), 
                         data=surv_loc_df, 
                         na.action = na.omit)


surv_reg_fixed=update(suv_model_simple, ~ . + Reg) ##just region as fixed effect
summary(surv_reg_fixed)

##################################
## INTERACTION OF F WITH REGION ##
##################################

surv_froh_inter=update(suv_model_simple, ~ . -FROH + (Reg*FROH)) #interaction between region and froh
summary(surv_froh_inter)

# FROH                    -18.652426   5.428689  -3.436 0.000591 ***
# RegLA:FROH                5.588950   7.767775   0.720 0.471830    
# RegMG:FROH                7.106856   7.116497   0.999 0.317967    
# RegNG:FROH                3.186312   6.875532   0.463 0.643058    
# RegSG:FROH               14.986234   8.374744   1.789 0.073541 .  
# RegSI:FROH                6.635056   7.199304   0.922 0.356724    

## emmtrends 
emtrends(surv_froh_inter, pairwise ~ Reg, var = "FROH")

# $emtrends
# Reg FROH.trend   SE  df asymp.LCL asymp.UCL
# IM      -18.65 5.43 Inf     -29.3     -8.01
# LA      -13.06 5.62 Inf     -24.1     -2.06
# MG      -11.55 4.65 Inf     -20.6     -2.44
# NG      -15.47 4.27 Inf     -23.8     -7.10
# SG       -3.67 6.38 Inf     -16.2      8.83
# SI      -12.02 4.78 Inf     -21.4     -2.65
# 
# Results are averaged over the levels of: Sex, MotherStatus 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast estimate   SE  df z.ratio p.value
# IM - LA    -5.589 7.77 Inf  -0.720  0.9796
# IM - MG    -7.107 7.12 Inf  -0.999  0.9185
# IM - NG    -3.186 6.88 Inf  -0.463  0.9973
# IM - SG   -14.986 8.37 Inf  -1.789  0.4725
# IM - SI    -6.635 7.20 Inf  -0.922  0.9411
# LA - MG    -1.518 7.26 Inf  -0.209  0.9999
# LA - NG     2.403 7.02 Inf   0.342  0.9994
# LA - SG    -9.397 8.49 Inf  -1.107  0.8785
# LA - SI    -1.046 7.36 Inf  -0.142  1.0000
# MG - NG     3.921 6.29 Inf   0.623  0.9894
# MG - SG    -7.879 7.88 Inf  -0.999  0.9183
# MG - SI     0.472 6.69 Inf   0.071  1.0000
# NG - SG   -11.800 7.68 Inf  -1.536  0.6407
# NG - SI    -3.449 6.43 Inf  -0.537  0.9947
# SG - SI     8.351 7.98 Inf   1.047  0.9021
# 
# Results are averaged over the levels of: Sex, MotherStatus 
# P value adjustment: tukey method for comparing a family of 6 estimates 

emmip(surv_froh_inter, Reg ~ FROH, cov.reduce = range,CIarg = list(lwd = 2, alpha = 0.5)) +
  theme_bw()


#plot predictions
inter=plot(predict_response(surv_froh_inter, terms = c("FROH[all]","Reg")),show.title=FALSE, line.size=1, colors="metro")+
  labs(x = expression(F["ROH"]), y = "Juvenile survival probability", colour = "Spatial \nregion")+
  theme(text = element_text(size = 15)) +
  xlim(0,0.2)
inter

inter_extended_surv=plot(predict_response(surv_froh_inter, terms = c("FROH[all]","Reg")),show.title=FALSE, line.size=1, colors="metro")+
  labs(x = expression(F["ROH"]), y = "Predicted juvenile survival probability", colour = "Spatial \nregion")+
  theme(text = element_text(size = 15)) 
inter_extended_surv




# only testing slope difference between 0-0.2 FROH because below this is where there are very few ids 

# test_predictions(surv_froh_inter,c("FROH[0,0.05,0.1,0.15,0.2]","Reg"), p_adjust = "bonferroni")

# (Average) Linear trend for FROH

# Reg   | Contrast |       95% CI |      p
# ----------------------------------------
# IM-LA |    -0.96 | -2.93,  1.01 | > .999
# IM-MG |    -1.60 | -3.30,  0.09 | 0.962 
# IM-NG |    -0.39 | -1.88,  1.10 | > .999
# IM-SG |    -2.72 | -4.30, -1.15 | 0.010 <<<<<
# IM-SI |    -0.82 | -2.57,  0.93 | > .999
# LA-MG |    -0.64 | -2.77,  1.49 | > .999
# LA-NG |     0.57 | -1.41,  2.56 | > .999
# LA-SG |    -1.76 | -3.81,  0.29 | > .999
# LA-SI |     0.14 | -2.06,  2.34 | > .999
# MG-NG |     1.22 | -0.50,  2.93 | > .999
# MG-SG |    -1.12 | -2.90,  0.66 | > .999
# MG-SI |     0.78 | -1.19,  2.75 | > .999
# NG-SG |    -2.33 | -3.93, -0.74 | 0.063 <<<<<
# NG-SI |    -0.44 | -2.22,  1.35 | > .999
# SG-SI |     1.90 |  0.05,  3.75 | 0.661 
# 
# Contrasts are presented as probabilities.


#test_predictions(surv_froh_inter,c("FROH[0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325]","Reg"), p_adjust = "bonferroni")

# (Average) Linear trend for FROH

# Reg   | Contrast |       95% CI |      p
# ----------------------------------------
#   IM-LA |    -0.40 | -1.87,  1.06 | > .999
# IM-MG |    -0.85 | -2.41,  0.72 | > .999
# IM-NG |    -0.15 | -0.98,  0.68 | > .999
# IM-SG |    -2.20 | -3.77, -0.63 | 0.090 
# IM-SI |    -0.51 | -1.64,  0.63 | > .999
# LA-MG |    -0.44 | -2.43,  1.54 | > .999
# LA-NG |     0.26 | -1.25,  1.76 | > .999
# LA-SG |    -1.80 | -3.79,  0.20 | > .999
# LA-SI |    -0.10 | -1.80,  1.60 | > .999
# MG-NG |     0.70 | -0.91,  2.30 | > .999
# MG-SG |    -1.35 | -3.42,  0.71 | > .999
# MG-SI |     0.34 | -1.46,  2.14 | > .999
# NG-SG |    -2.05 | -3.66, -0.44 | 0.187 
# NG-SI |    -0.36 | -1.56,  0.85 | > .999
# SG-SI |     1.69 | -0.10,  3.49 | 0.963 

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

save(inter_extended_surv, file = "Deer_spatial_variation_ID/plots/surv_IDxE_plots.RData")

