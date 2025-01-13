######################
## Annual Denisty ####
######################

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

surv_loc_df=read.table("Deer_spatial_variation_ID/survival_loc_2024.txt", sep = ",", header = TRUE)%>%
  select(Code, MumCode, BirthYear,Sex, N, E, Reg)%>%
  filter(Sex!=3)%>%
  na.omit()

density=read.csv("Deer_spatial_variation_ID/Density_GrazeType.csv", header=T, stringsAsFactors = F)%>%
  rename(MumCode=Name, BirthYear=Year)

loc_den_graze=surv_loc_df%>%
  inner_join(density)%>%
  unique()%>%
  na.omit()

den_simple=glmmTMB(AnnualDensity~ Sex + Reg, 
                   family="gaussian", 
                   data=loc_den_graze, 
                   na.action = na.omit,
)


summary(den_simple)

den_pred_plot=
  ggplot()+
  theme_bw()+
  scale_color_manual(values = c("#f0a30a" ,"#a20025","#00aba9","chocolate1", "#60a917","#647687"))+
  labs(x="Spatial region", y="Derived population denisty metric")+
  theme(text = element_text(size = 18),legend.position = "none")+
  geom_boxplot(data=loc_den_graze, aes(Reg, AnnualDensity, colour=Reg, alpha=0.1),
               inherit.aes = F, position=position_nudge(x=0.2), width=0.1)


den_pred_plot


graze_simple=glmmTMB(GrazeType~ Sex + Reg, 
                   family="gaussian", 
                   data=loc_den_graze, 
                   na.action = na.omit,
)


summary(graze_simple)

graz_pred_plot=
  ggplot()+
  theme_bw()+
  scale_color_manual(values = c("#f0a30a" ,"#a20025","#00aba9","chocolate1", "#60a917","#647687"))+
  labs(x="Spatial region", y="Grazing quality")+
  theme(text = element_text(size = 18),legend.position = "none")+
  geom_boxplot(data=loc_den_graze, aes(Reg, GrazeType, colour=Reg, alpha=0.1),
               inherit.aes = F, position=position_nudge(x=0.2), width=0.1)


graz_pred_plot

box=graz_pred_plot+den_pred_plot+plot_annotation(tag_levels = 'A')

ggsave(box,
       file = "Deer_spatial_variation_ID/plots/SUPP_boxplots_graze_popDen.png",
       width = 7,
       height = 6, 
       bg = "white"
)


###################
### using INLA ###
##################

loc_den_graze_2=loc_den_graze%>%mutate(AnnualDensity=AnnualDensity*1000) #rememnber is in e-03 now


IM1_den  <- inla(AnnualDensity~ Sex + Reg, 
                 family = "gaussian",
                 data = loc_den_graze_2,
                 control.compute = list(dic=TRUE)) 


summary(IM1_den)


rum_outline=read.csv("PhD/Chapter_5_spatial_ID_x_E/Spatial_var_inbreeding/INLA/RumBoundary.csv")%>%
  rename(E = Easting, N = Northing) 


N=nrow(rum_outline)
rum_line_rev=rum_outline[N:1, c("E","N")]

Mesh=inla.mesh.2d(loc.domain= rum_outline, 
                  max.edge=2, #probs use 1 for actual model
                  boundary=
                    inla.mesh.segment(rum_line_rev))

plot(Mesh, asp=1)

Locations=cbind(loc_den_graze$E, loc_den_graze$N)#locations of ids

A=inla.spde.make.A(Mesh, loc=Locations)

#define SPDE 
spde=inla.spde2.matern(Mesh, alpha = 2) # would need to adjust for time series data

#define spatial field
w.index=inla.spde.make.index(name = 'w', n.spde=spde$n.spde, n.group = 1, n.repl = 1)

#make model matrix
N <- nrow(loc_den_graze)
X0=data.frame(Intercept = rep(1, N),
              Sex = loc_den_graze$Sex, 
              Reg= loc_den_graze$Reg)

x=as.data.frame(X0)


stackfit_den=inla.stack(
  data=list(y=loc_den_graze_2$AnnualDensity), 
  A = list(1,  A), 
  effects=list(
    X=x,
    w=w.index)
)

#define SPDE 


IM_spde_den  <- inla(y~ -1 + Intercept+Sex + Reg+ f(w, model=spde), 
                     family = "gaussian",
                     data=inla.stack.data(stackfit_den), 
                     control.compute = list(dic=TRUE),
                     control.predictor = list(
                       A=inla.stack.A(stackfit_den))) 


summary(IM_spde_den)

inla_den_plot=ggField(IM_spde_den, Mesh)+
  labs(fill = "Derived \npopulation \ndensity")+
  theme_bw()+
  scale_fill_discrete_sequential(labels=c("Low", " "," "," "," "," "," "," ","High"), 
                                 palette = "Purples", rev=T)+
  theme(text = element_text(size = 18),
        legend.title=element_text(size=rel(0.8))) +
  annotate("segment", x = 1372, xend = 1382, y = 8000, yend = 8000, colour = "black", linewidth = 1) +
  annotate("text" ,x = 1377, y = 8001.5, label = "1km")


inla_den_plot

SpatialList <- list(IM1_den,  IM_spde_den)
sapply(SpatialList, function(f) f$dic$dic)
INLADICFig(SpatialList)+theme_classic()


### and for grazing quality too ###

stackfit2=inla.stack(
  data=list(y=loc_den_graze$GrazeType), 
  A = list(1,  A), 
  effects=list(
    X=x,
    w=w.index)
)

IM_spde_graze  <- inla(y~ -1 + Intercept+Sex + Reg+ f(w, model=spde), 
                       family = "gaussian",
                       data=inla.stack.data(stackfit2), 
                       control.compute = list(dic=TRUE),
                       control.predictor = list(
                         A=inla.stack.A(stackfit2))) 


summary(IM_spde_graze)

inla_graze_plot=ggField(IM_spde_graze, Mesh)+
  labs(fill = "Grazing Quality")+
  theme_bw()+
  scale_fill_discrete_sequential(palette = "Greens", rev=T, 
                                 labels=c("Low", " "," "," "," "," "," "," ","High"))+
  theme(text = element_text(size = 18),
        legend.title=element_text(size=rel(0.8))) +
  annotate("segment", x = 1372, xend = 1382, y = 8000, yend = 8000, colour = "black", linewidth = 1) +
  annotate("text" ,x = 1377, y = 8001.5, label = "1km")


inla_graze_plot

IM1  <- inla(GrazeType~ Sex + Reg, 
             
             family = "gaussian",
             data = loc_den_graze,
             control.compute = list(dic=TRUE)) 

SpatialList <- list(IM1,  IM_spde_graze)
sapply(SpatialList, function(f) f$dic$dic)
INLADICFig(SpatialList)+theme_classic()




inla=inla_graze_plot+inla_den_plot+plot_annotation(tag_levels = 'A')
inla

ggsave(inla,
       file = "Deer_spatial_variation_ID/plots/SUPP_INLA_graze_popDen.png",
       width = 10,
       height = 5, 
       bg = "white"
)
