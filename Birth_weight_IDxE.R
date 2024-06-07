library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(patchwork)
library(emmeans)
library(INLA)
library(inlabru)
library(ggregplot)

setwd("/Volumes/Seagate Por")

surv_loc_df=read.table("Deer_spatial_variation_ID/survival_loc.txt", sep = ",", header = TRUE)%>%
  select(-BirthWt)

birth_wt<-read.csv("Deer_spatial_variation_ID/sys_BirthWt.csv") 

bw_df=surv_loc_df%>%full_join(birth_wt)%>%na.omit
head(bw_df)
table(bw_df$Reg)

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


inter=plot(predict_response(bw_reg_inter, terms = c("FROH[all]","Reg")),show.title=FALSE, line.size=1, colors="metro")+
  labs(x = expression(F["ROH"]), y = "Birth Weight (kg)", colour = "Spatial \nregion")+
  theme(text = element_text(size = 15)) 
 # xlim(0,0.2)+
inter



# only testing slope difference between 0-0.2 FROH because below this is where there are very few ids 
test_predictions(bw_reg_inter,c("FROH[0,0.05,0.1,0.15,0.2]","Reg"), p_adjust = "bonferroni")

# Reg   | Contrast |       95% CI |      p
# ----------------------------------------
# IM-LA |    -3.98 |  -9.75, 1.79 | > .999
# IM-MG |    -2.17 |  -7.47, 3.13 | > .999
# IM-NG |    -2.22 |  -6.67, 2.24 | > .999
# IM-SG |    -5.67 | -11.98, 0.64 | > .999
# IM-SI |    -4.40 |  -9.24, 0.44 | > .999
# LA-MG |     1.81 |  -4.18, 7.80 | > .999
# LA-NG |     1.77 |  -3.49, 7.02 | > .999
# LA-SG |    -1.68 |  -8.59, 5.22 | > .999
# LA-SI |    -0.42 |  -6.04, 5.21 | > .999
# MG-NG |    -0.04 |  -4.82, 4.73 | > .999
# MG-SG |    -3.49 | -10.01, 3.02 | > .999
# MG-SI |    -2.23 |  -7.39, 2.93 | > .999
# NG-SG |    -3.45 |  -9.32, 2.42 | > .999
# NG-SI |    -2.18 |  -6.47, 2.11 | > .999
# SG-SI |     1.27 |  -4.91, 7.44 | > .999
# 
