library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(patchwork)
library(emmeans)
library(INLA)
library(inlabru)
library(ggregplot)

setwd("/Volumes/Seagate_HD/")

## base models + region (but no FROH)
# to show that the effect is independent from inbreeding

bw_df=read.table("Deer_spatial_variation_ID/birthweight_loc_2024_new_calves.txt", sep = ",", header = TRUE)

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

summary(bw_model_simple)

bw_reg_fixed=update(bw_model_simple, ~ . + Reg) ##just region as fixed effect
summary(bw_reg_fixed)


## SURVIVAL ##
surv_loc_df=read.table("Deer_spatial_variation_ID/survival_loc_2024_new_calves.txt", sep = ",", header = TRUE)%>%
  select(-MumFROH, -BirthWt)%>%
  filter(!E>1385)%>%
  filter(!N<7997.5) %>%#removing ids with no known region or ~10 ids with outside the limits of study area
  na.omit()


head(surv_loc_df)

surv_loc_df$Code=as.factor(surv_loc_df$Code)
surv_loc_df$MumCode=as.factor(surv_loc_df$MumCode)
surv_loc_df$BirthYear=as.factor(surv_loc_df$BirthYear)
surv_loc_df$Sex=as.factor(surv_loc_df$Sex)
surv_loc_df$MotherStatus=as.factor(surv_loc_df$MotherStatus)
surv_loc_df$Reg=as.factor(surv_loc_df$Reg)
## base model of juvenile survival
suv_model_simple=glmmTMB(juvenile_survival~ 1+ Sex + MotherStatus + mum_age+mum_age_sq+Day_seq+
                           (1|BirthYear)+(1|MumCode), 
                         family=binomial, 
                         data=surv_loc_df, 
                         na.action = na.omit)

summary(suv_model_simple)


surv_reg_fixed=update(suv_model_simple, ~ . + Reg) ##just region as fixed effect
summary(surv_reg_fixed)
