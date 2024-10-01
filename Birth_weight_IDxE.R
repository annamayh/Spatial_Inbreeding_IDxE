library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(patchwork)
library(emmeans)


setwd("/Volumes/Seagate_HD/")

bw_df=read.table("Deer_spatial_variation_ID/birthweight_loc_2024.txt", sep = ",", header = TRUE)

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

bw_reg_inter=update(bw_model_simple, ~ . + Reg*FROH) ##just region as fixed effect
summary(bw_reg_inter)

emtrends(bw_reg_inter, pairwise ~ Reg, var = "FROH")


# $emtrends
# Reg FROH.trend   SE   df lower.CL upper.CL
# IM      -4.892 1.67 2362    -8.17    -1.62
# LA      -3.079 2.15 2362    -7.29     1.13
# MG      -4.590 1.78 2362    -8.08    -1.10
# NG      -4.044 1.33 2362    -6.66    -1.43
# SG      -0.139 2.34 2362    -4.74     4.46
# SI      -1.863 1.64 2362    -5.08     1.35
# 
# Results are averaged over the levels of: Sex, MotherStatus 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast estimate   SE   df t.ratio p.value
# IM - LA    -1.813 2.71 2362  -0.668  0.9854
# IM - MG    -0.302 2.44 2362  -0.124  1.0000
# IM - NG    -0.849 2.13 2362  -0.398  0.9987
# IM - SG    -4.753 2.88 2362  -1.653  0.5636
# IM - SI    -3.029 2.33 2362  -1.298  0.7864
# LA - MG     1.511 2.78 2362   0.543  0.9944
# LA - NG     0.965 2.52 2362   0.383  0.9989
# LA - SG    -2.940 3.18 2362  -0.926  0.9399
# LA - SI    -1.216 2.70 2362  -0.450  0.9977
# MG - NG    -0.546 2.23 2362  -0.245  0.9999
# MG - SG    -4.451 2.94 2362  -1.514  0.6555
# MG - SI    -2.727 2.43 2362  -1.122  0.8722
# NG - SG    -3.905 2.69 2362  -1.450  0.6965
# NG - SI    -2.180 2.13 2362  -1.025  0.9097
# SG - SI     1.724 2.86 2362   0.602  0.9909
# 
# Results are averaged over the levels of: Sex, MotherStatus 
# P value adjustment: tukey method for comparing a family of 6 estimates

emmip(bw_reg_inter, Reg ~ FROH, cov.reduce = range,CIarg = list(lwd = 2, alpha = 0.5)) +
  theme_bw()

## predictions from model

inter=plot(predict_response(bw_reg_inter, terms = c("FROH[all]","Reg")),show.title=FALSE, line.size=1, colors="metro")+
  labs(x = expression(F["ROH"]), y = "Birth Weight (kg)", colour = "Spatial \nregion")+
  theme(text = element_text(size = 15)) +
  xlim(0,0.2)

inter

inter_extended_bw=plot(predict_response(bw_reg_inter, terms = c("FROH[all]","Reg")),show.title=FALSE, line.size=1, colors="metro")+
  labs(x = expression(F["ROH"]), y = "Predicted birth weight (kg)", colour = "Spatial \nregion")+
  theme(text = element_text(size = 15)) 
inter_extended_bw


# only testing slope difference between 0-0.2 FROH because below this is where there are very few ids 
#test_predictions(bw_reg_inter,c("FROH[0,0.05,0.1,0.15,0.2]","Reg"), p_adjust = "bonferroni")

# (Average) Linear trend for FROH
# 
# Reg   | Contrast |       95% CI |      p
# ----------------------------------------
#   IM-LA |    -1.81 |  -7.13, 3.51 | > .999
# IM-MG |    -0.30 |  -5.08, 4.47 | > .999
# IM-NG |    -0.85 |  -5.02, 3.33 | > .999
# IM-SG |    -4.75 | -10.39, 0.89 | > .999
# IM-SI |    -3.03 |  -7.60, 1.55 | > .999
# LA-MG |     1.51 |  -3.94, 6.96 | > .999
# LA-NG |     0.96 |  -3.97, 5.90 | > .999
# LA-SG |    -2.94 |  -9.16, 3.28 | > .999
# LA-SI |    -1.22 |  -6.51, 4.08 | > .999
# MG-NG |    -0.55 |  -4.92, 3.82 | > .999
# MG-SG |    -4.45 | -10.21, 1.31 | > .999
# MG-SI |    -2.73 |  -7.49, 2.03 | > .999
# NG-SG |    -3.90 |  -9.19, 1.38 | > .999
# NG-SI |    -2.18 |  -6.35, 1.99 | > .999
# SG-SI |     1.72 |  -3.89, 7.34 | > .999
# 



emtrends(bw_reg_inter, pairwise ~ Reg, var = "FROH")

emmip(bw_reg_inter, Reg ~ FROH, cov.reduce = range)



ggsave(inter,
       file = "Deer_spatial_variation_ID/plots/IDxE_birthWt.png",
       width = 7,
       height = 6, 
       bg = "white"
)


ggsave(inter_extended,
       file = "Deer_spatial_variation_ID/plots/IDxE_birthWt_SUPP.png",
       width = 7,
       height = 6, 
       bg = "white"
)



save(inter_extended_bw, file = "Deer_spatial_variation_ID/plots/bw_IDxE.RData")

