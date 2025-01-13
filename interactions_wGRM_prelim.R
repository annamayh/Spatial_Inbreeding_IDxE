library(tidyverse)
library(brms)
library(corpcor)

setwd("/Volumes/Seagate_HD/")

surv_loc_df=read.table("Deer_spatial_variation_ID/survival_loc_2024.txt", sep = ",", header = TRUE)#%>%

grm=readRDS(file = "PhD/Deer_SNP_data/GRM_matrix_names.RDS")
grm_filt=grm[rownames(grm)%in%surv_loc_df[["Code"]],colnames(grm)%in%surv_loc_df[["Code"]]]##filtering grm for ids only in df
grm_filt_pd <- make.positive.definite(grm_filt)
GRM <- as(grm_filt_pd, "dgCMatrix")

surv_loc_df$Code=as.factor(surv_loc_df$Code)
surv_loc_df$MumCode=as.factor(surv_loc_df$MumCode)
surv_loc_df$BirthYear=as.factor(surv_loc_df$BirthYear)
surv_loc_df$Sex=as.factor(surv_loc_df$Sex)
surv_loc_df$MotherStatus=as.factor(surv_loc_df$MotherStatus)
surv_loc_df$Reg=as.factor(surv_loc_df$Reg)

#standard model with no interactions and no GRM
suv_model_simple=brm(juvenile_survival~ 1+ Sex + MotherStatus + mum_age+mum_age_sq+Day_seq+FROH+
                     (1|BirthYear)+(1|MumCode), 
                     family=bernoulli, 
                     data=surv_loc_df, 
                     chains=4, 
                     cores = 4
)

summary(suv_model_simple)

## default number of itts actually does well
suv_model_wGRM=brm(juvenile_survival~ 1+ Sex + MotherStatus + mum_age+mum_age_sq+Day_seq+FROH+
                      (1|gr(Code, cov=Amat)) + (1|BirthYear)+(1|MumCode), 
                      family=bernoulli, 
                      data=surv_loc_df, 
                      data2 = list(Amat = GRM),
                     chains=4, 
                     cores = 4
                     )



suv_model_plus_inter=brm(juvenile_survival~ 1+ Sex + MotherStatus + mum_age+mum_age_sq+Day_seq+(FROH*Reg)+
                       (1|BirthYear)+(1|MumCode), 
                     family=bernoulli, 
                     data=surv_loc_df, 
                     chains=4, 
                     cores = 4, 
                     iter = 6000, 
                     warmup = 2000, 
                     thin = 5
)
summary(suv_model_plus_inter)


suv_model_plus_inter_andGRM=brm(juvenile_survival~ 1+ Sex + MotherStatus + mum_age+mum_age_sq+Day_seq+(FROH*Reg)+
                          (1|gr(Code, cov=Amat))+(1|BirthYear)+(1|MumCode), 
                         family=bernoulli, 
                         data=surv_loc_df, 
                         data2 = list(Amat = GRM),
                         chains=4, 
                         cores = 4
)
summary(suv_model_plus_inter_andGRM)

plot(suv_model_simple)


newdata <- data.frame(#Reg = factor(c("NG")),
                      FROH = c(0,0.05, 0.1))

pred=brms::predict(suv_model_plus_inter_andGRM, newdata = newdata)
pred
