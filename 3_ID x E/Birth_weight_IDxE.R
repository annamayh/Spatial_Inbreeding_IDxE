library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(patchwork)
library(emmeans)

#updated 25/11/24

setwd("/Volumes/Seagate_HD/")

bw_df=read.table("Deer_spatial_variation_ID/birthweight_loc_2024_new_calves.txt", sep = ",", header = TRUE)%>%
  mutate(FROH_trans=sqrt(FROH))


bw_df$Code=as.factor(bw_df$Code)
bw_df$MumCode=as.factor(bw_df$MumCode)
bw_df$BirthYear=as.factor(bw_df$BirthYear)
bw_df$Sex=as.factor(bw_df$Sex)
bw_df$Reg=as.factor(bw_df$Reg)

## Region as categorical ##

bw_model_simple=glmmTMB(CaptureWt~ Sex + AgeHrs+ MotherStatus+mum_age+mum_age_sq+Day_seq+
                          (1|BirthYear)+ (1|MumCode), 
                        family=gaussian(), 
                        data=bw_df, 
                        na.action = na.omit,
)

bw_reg_inter=update(bw_model_simple, ~ . + Reg*FROH_trans) ##
summary(bw_reg_inter)

emtrends(bw_reg_inter, pairwise ~ Reg, var = "FROH_trans")
# $emtrends
# Reg FROH_trans.trend    SE   df lower.CL upper.CL
# IM            -3.517 1.020 2478    -5.52  -1.5102
# LA            -2.237 1.120 2478    -4.44  -0.0339
# MG            -1.738 0.897 2478    -3.50   0.0216
# NG            -2.254 0.703 2478    -3.63  -0.8746
# SG            -0.268 1.170 2478    -2.57   2.0306
# SI            -1.185 0.873 2478    -2.90   0.5259
# 
# Results are averaged over the levels of: Sex, MotherStatus 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast estimate   SE   df t.ratio p.value
# IM - LA   -1.2800 1.52 2478  -0.844  0.9592
# IM - MG   -1.7789 1.36 2478  -1.311  0.7793
# IM - NG   -1.2627 1.24 2478  -1.019  0.9117
# IM - SG   -3.2485 1.55 2478  -2.091  0.2921
# IM - SI   -2.3315 1.34 2478  -1.735  0.5083
# LA - MG   -0.4989 1.44 2478  -0.348  0.9993
# LA - NG    0.0173 1.32 2478   0.013  1.0000
# LA - SG   -1.9685 1.62 2478  -1.213  0.8309
# LA - SI   -1.0515 1.42 2478  -0.740  0.9769
# MG - NG    0.5162 1.14 2478   0.454  0.9976
# MG - SG   -1.4696 1.48 2478  -0.990  0.9213
# MG - SI   -0.5526 1.26 2478  -0.440  0.9979
# NG - SG   -1.9858 1.37 2478  -1.455  0.6933
# NG - SI   -1.0688 1.13 2478  -0.948  0.9338
# SG - SI    0.9170 1.46 2478   0.628  0.9890
# Results are averaged over the levels of: Sex, MotherStatus 
# P value adjustment: tukey method for comparing a family of 6 estimates

## predictions from model



pred_bw=predict_response(bw_reg_inter, terms = c("FROH_trans[all]","Reg"))%>%
  as.data.frame()%>%
  mutate(FROH=x^2)## tranforming back into FROH

cols = c("SI"="#f0a30a" ,"IM" = "#a20025", "LA" = "#1ba1e2", "NG" = "#fa6800", "MG" = "#008a00", "SG" = "#647687")

inter_bw=ggplot(pred_bw, aes(x=FROH, y=predicted, ymin=conf.low, ymax=conf.high,
                                 group = group, colour = group, fill=group))+
  geom_line(linewidth=1)+
  geom_ribbon(alpha=0.15, colour = NA, show.legend = F)+
  labs(x = expression(F["ROH"]), y = "Predicted birth weight (kg)", colour = "Spatial \nregion")+
  theme_bw()+
  theme(text = element_text(size = 15)) +
  scale_color_manual(values = cols, breaks=c("SI","IM","LA", "NG", "MG", "SG"))+
  scale_fill_manual(values = cols)

inter_bw






# ggsave(inter,
#        file = "Deer_spatial_variation_ID/plots/IDxE_birthWt.png",
#        width = 7,
#        height = 6, 
#        bg = "white"
# )
# 
# 
# ggsave(inter_extended,
#        file = "Deer_spatial_variation_ID/plots/IDxE_birthWt_SUPP.png",
#        width = 7,
#        height = 6, 
#        bg = "white"
# )



save(inter_bw, file = "Deer_spatial_variation_ID/plots/bw_IDxE.RData")
# updated 25.11.23
