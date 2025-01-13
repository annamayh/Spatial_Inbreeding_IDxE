library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(patchwork)
library(emmeans)
library(INLA)
library(ggregplot)
library(inlabru)
library(colorspace)

setwd("/Volumes/Seagate_HD/")

surv_loc_df=read.table("Deer_spatial_variation_ID/survival_loc_2024.txt", sep = ",", header = TRUE)

density=read.csv("Deer_spatial_variation_ID/Density_GrazeType.csv", header=T, stringsAsFactors = F)%>%
  rename(MumCode=Name, BirthYear=Year)

loc_den_graze=surv_loc_df%>%
  inner_join(density)%>%
  unique()


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
                         data=loc_den_graze, 
                         na.action = na.omit)


surv_reg_fixed=update(suv_model_simple, ~ . + AnnualDensity) ##just region as fixed effect
summary(surv_reg_fixed)

surv_grazing=update(suv_model_simple, ~ . + GrazeType) ##just region as fixed effect
summary(surv_grazing)

##################################
## INTERACTION OF F WITH REGION ##
##################################

surv_froh_interDen=update(suv_model_simple, ~ . -FROH + (AnnualDensity*FROH)) #interaction between region and froh
summary(surv_froh_interDen)

surv_froh_inter=update(suv_model_simple, ~ . -FROH + (GrazeType*FROH)) #interaction between region and froh
summary(surv_froh_inter)
