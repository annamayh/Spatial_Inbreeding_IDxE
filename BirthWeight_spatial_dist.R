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

bw_reg_fixed=update(bw_model_simple, ~ . + Reg+FROH) ##just region as fixed effect
summary(bw_reg_fixed)


reg_pred_f=ggpredict(bw_reg_fixed, terms = c("Reg","Sex[1]","AgeHrs[0]"))%>%
  mutate(x = fct_relevel(x, "SI", "IM", "LA", "NG","MG",  "SG"))

reg_pred_m=ggpredict(bw_reg_fixed, terms = c("Reg","Sex[2]","AgeHrs[0]"))%>%
  mutate(x = fct_relevel(x, "SI", "IM", "LA", "NG","MG",  "SG"))

reg_pred=rbind(reg_pred_f,reg_pred_m)%>%
  arrange(x)%>%na.omit()

reg_pred

bw_reg=reg_pred%>%
  ggplot(aes(x=x, y=predicted, color=x, ymin=conf.low, ymax=conf.high, group=group))+
  geom_pointrange(linewidth=1, position = position_dodge(width=0.5))+
  theme_bw()+
  scale_color_manual(values = c("#f0a30a" ,"#a20025","#00aba9","chocolate1", "#60a917","#647687"))+
  labs(x="Spatial region", y="Predicted birth weight (kg)")+
  theme(text = element_text(size = 18),legend.position = "none")
bw_reg



### Region as continuous in INLA ###

IM2  <- inla(CaptureWt~ Sex + AgeHrs+ MotherStatus+FROH+
               f(BirthYear, model = 'iid')+f(MumCode, model = 'iid'), 
             
             family = "gaussian",
             data = bw_df,
             control.compute = list(dic=TRUE)) 


summary(IM2)


## add spatial field
rum_outline=read.csv("PhD/Chapter_5_spatial_ID_x_E/Spatial_var_inbreeding/INLA/RumBoundary.csv")%>%
  rename(E = Easting, N = Northing) 


N=nrow(rum_outline)
rum_line_rev=rum_outline[N:1, c("E","N")]

Mesh=inla.mesh.2d(loc.domain= rum_outline, 
                  max.edge=2, #probs use 1 for actual model
                  boundary=
                    inla.mesh.segment(rum_line_rev))

plot(Mesh, asp=1)

Locations=cbind(bw_df$E, bw_df$N)#locations of ids

A=inla.spde.make.A(Mesh, loc=Locations)

#define SPDE 
spde=inla.spde2.matern(Mesh, alpha = 2) # would need to adjust for time series data

#define spatial field
w.index=inla.spde.make.index(name = 'w', n.spde=spde$n.spde, n.group = 1, n.repl = 1)

#make model matrix
N <- nrow(bw_df)
X0=data.frame(Intercept = rep(1, N),
              Sex = bw_df$Sex, 
              MotherStatus=bw_df$MotherStatus, 
              AgeHrs=bw_df$AgeHrs, 
              FROH=bw_df$FROH,
              mum_age=bw_df$mum_age,
              mum_age_sq=bw_df$mum_age_sq,
              Day_seq=bw_df$Day_seq)

x=as.data.frame(X0)


## now fitting spde and birth year as random effects 
# have to re-do the stack
stackfit2=inla.stack(
  data=list(y=bw_df$CaptureWt), 
  A = list(1, 1, 1, A), 
  effects=list(
    X=x,
    BirthYear = bw_df$BirthYear,
    MumCode = bw_df$MumCode, # insert vectors of any random effects
    w=w.index)
)

#define SPDE 


IM_spde  <- inla(y~ -1 + Intercept+Sex + MotherStatus + AgeHrs+FROH+mum_age+mum_age_sq+Day_seq+
                   f(BirthYear, model = 'iid')+f(MumCode, model = 'iid')+ f(w, model=spde), 
                 family = "gaussian",
                 data=inla.stack.data(stackfit2), 
                 control.compute = list(dic=TRUE),
                 control.predictor = list(
                   A=inla.stack.A(stackfit2))) 


summary(IM_spde)


List <- list(IM2, IM_spde)
sapply(List, function(f) f$dic$dic)
INLADICFig(List)

library(colorspace)

inla_bw_plot=ggField(IM_spde, Mesh)+
  labs(fill = "Birth weight (kg) \n(as deviation \nfrom mean)")+
  theme_bw()+
  scale_fill_discrete_sequential(palette = "Burg", rev=FALSE)+
  theme(text = element_text(size = 18),
        legend.title=element_text(size=rel(0.7)))


inla_bw_plot


birth_weight=bw_reg+inla_bw_plot+plot_annotation(tag_levels = 'A')

save(bw_reg,inla_bw_plot, file = "Deer_spatial_variation_ID/plots/bw_plots.RData")


# ggsave(birth_weight,
#        file = "Deer_spatial_variation_ID/plots/bw_spatial_var.png",
#        width = 11,
#        height = 10,
#        dpi=1000)
