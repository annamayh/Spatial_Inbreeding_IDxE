library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(patchwork)
library(emmeans)
library(INLA)
library(inlabru)
library(ggregplot)


setwd("/Volumes/Seagate_HD/")

surv_loc_df=read.table("Deer_spatial_variation_ID/survival_loc_2024.txt", sep = ",", header = TRUE)%>%
  select(-BirthWt)%>%
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


IM1  <- inla(juvenile_survival~1+Sex + MotherStatus + mum_age+mum_age_sq+Day_seq+FROH+
               f(BirthYear, model = 'iid')+f(MumCode, model = 'iid'), 
             family = "binomial",
             data = surv_loc_df,
             control.compute = list(dic=TRUE)) 


summary(IM1)



IM2  <- inla(juvenile_survival~1+Sex + MotherStatus + mum_age+mum_age_sq+Day_seq+FROH+MumFROH+
               f(BirthYear, model = 'iid')+f(MumCode, model = 'iid'), 
             family = "binomial",
             data = surv_loc_df,
             control.compute = list(dic=TRUE)) 


summary(IM2)

##### now add in spatial variation #####
rum_outline=read.csv("PhD/Chapter_5_spatial_ID_x_E/Spatial_var_inbreeding/INLA/RumBoundary.csv")%>%
  rename(E = Easting, N = Northing) 

N=nrow(rum_outline)
rum_line_rev=rum_outline[N:1, c("E","N")]

Mesh=inla.mesh.2d(loc.domain= rum_outline, 
                  max.edge=2, #probs use 1 for actual model
                  boundary=
                    inla.mesh.segment(rum_line_rev))



Locations=cbind(surv_loc_df$E, surv_loc_df$N)#locations of ids
#make A mtrix
A=inla.spde.make.A(Mesh, loc=Locations)
dim(A)

#define SPDE 
spde=inla.spde2.matern(Mesh, alpha = 2) # would need to adjust for time series data

#define spatial field
w.index=inla.spde.make.index(name = 'w', n.spde=spde$n.spde, n.group = 1, n.repl = 1)

#make model matrix
N <- nrow(surv_loc_df)
X0=data.frame(Intercept = rep(1, N),
              FROH = surv_loc_df$FROH, 
              Sex = surv_loc_df$Sex, 
              MotherStatus = surv_loc_df$MotherStatus, 
              mum_age =surv_loc_df$mum_age,
              mum_age_sq = surv_loc_df$mum_age_sq, 
              Day_seq=surv_loc_df$Day_seq)
x=as.data.frame(X0)

#make stack
stackfit=inla.stack(
  tag="Fit",
  data=list(y=surv_loc_df$juvenile_survival), 
  A = list(A, 1, 1, 1), #for the 6 fixed and random effects I have 
  effects=list(
    w=w.index,
    X=x,
    BirthYear = surv_loc_df$BirthYear, 
    MumCode = surv_loc_df$MumCode)
)


IM1_spde  <- inla(y~ -1 + Intercept+Sex + MotherStatus + mum_age+mum_age_sq+Day_seq+FROH+
                   f(BirthYear, model = 'iid')+f(MumCode, model = 'iid')+ f(w, model=spde), 
                 family = "binomial",
                 data=inla.stack.data(stackfit), 
                 control.compute = list(dic=TRUE),
                 control.predictor = list(
                   A=inla.stack.A(stackfit))) 


summary(IM1_spde)


## re-make stack etc.

X0_2=data.frame(Intercept = rep(1, N),
              FROH = surv_loc_df$FROH, 
              Sex = surv_loc_df$Sex, 
              MotherStatus = surv_loc_df$MotherStatus, 
              mum_age =surv_loc_df$mum_age,
              mum_age_sq = surv_loc_df$mum_age_sq, 
              Day_seq=surv_loc_df$Day_seq,
              MumFROH = surv_loc_df$MumFROH)
x2=as.data.frame(X0_2)

#make stack
stackfit2=inla.stack(
  tag="Fit",
  data=list(y=surv_loc_df$juvenile_survival), 
  A = list(A, 1, 1, 1), #for the 6 fixed and random effects I have 
  effects=list(
    w=w.index,
    X=x2,
    BirthYear = surv_loc_df$BirthYear, 
    MumCode = surv_loc_df$MumCode)
)


IM2_spde  <- inla(y~ -1 + Intercept+Sex + MotherStatus + mum_age+mum_age_sq+Day_seq+FROH+MumFROH+
                    f(BirthYear, model = 'iid')+f(MumCode, model = 'iid')+ f(w, model=spde), 
                  family = "binomial",
                  data=inla.stack.data(stackfit2), 
                  control.compute = list(dic=TRUE),
                  control.predictor = list(
                  A=inla.stack.A(stackfit2))) 


summary(IM2_spde)





SpatialList <- list(IM1, IM2, IM1_spde, IM2_spde)
sapply(SpatialList, function(f) f$dic$dic)
INLADICFig(SpatialList, ModelNames = c('IM1', 'IM2', 'IM1_spde', 'IM2_spde'))+theme_classic()


inla_surv_plot=ggField(IM_spde, Mesh)+
  labs(fill = "Juvenile survival \n(as untransformed \ndevaition from mean)")+
  theme_bw()+
  scale_fill_discrete_sequential(palette = "Oranges", rev=FALSE)+
  theme(text = element_text(size = 18),
        legend.title=element_text(size=rel(0.7))) +
  annotate("segment", x = 1372, xend = 1382, y = 8000, yend = 8000, colour = "black", linewidth = 1) +
  annotate("text" ,x = 1377, y = 8001.5, label = "1km")

inla_surv_plot



surv_reg+inla_surv_plot+plot_annotation(tag_levels = 'A')
save(surv_reg,inla_surv_plot, file = "Deer_spatial_variation_ID/plots/surv_plots.RData")
