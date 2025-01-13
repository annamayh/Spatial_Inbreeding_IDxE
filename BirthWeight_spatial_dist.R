library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(patchwork)
library(emmeans)
library(INLA)
library(inlabru)
library(ggregplot)
library(effsize)

setwd("/Volumes/Seagate_HD/")

#upadted 25.11.24
bw_df=read.table("Deer_spatial_variation_ID/birthweight_loc_2024_new_calves.txt", sep = ",", header = TRUE)

mean(bw_df$CaptureWt)
sd(bw_df$CaptureWt)

bw_df%>%
  group_by(Reg)%>%
  summarise(mean=mean(CaptureWt))


bw_df$Code=as.factor(bw_df$Code)
bw_df$MumCode=as.factor(bw_df$MumCode)
bw_df$BirthYear=as.factor(bw_df$BirthYear)
bw_df$Sex=as.factor(bw_df$Sex)
bw_df$Reg=as.factor(bw_df$Reg)

cohen.d(bw_df$CaptureWt[bw_df$Reg == "SI"], bw_df$CaptureWt[bw_df$Reg == "SG"])
cohen.d(bw_df$CaptureWt[bw_df$Reg == "SI"], bw_df$CaptureWt[bw_df$Reg == "LA"])
cohen.d(bw_df$CaptureWt[bw_df$Reg == "SI"], bw_df$CaptureWt[bw_df$Reg == "NG"])
cohen.d(bw_df$CaptureWt[bw_df$Reg == "SI"], bw_df$CaptureWt[bw_df$Reg == "MG"])


## Region as categorical ##

bw_model_simple=glmmTMB(CaptureWt~ Sex + AgeHrs+ MotherStatus+mum_age+mum_age_sq+Day_seq+FROH+
                          (1|BirthYear)+ (1|MumCode), 
                        family=gaussian(), 
                        data=bw_df, 
                        na.action = na.omit,
)

bw_reg_fixed=update(bw_model_simple, ~ . + Reg) ##just region as fixed effect
summary(bw_reg_fixed)

anova(bw_model_simple,bw_reg_fixed) ## region fixed is a better fit i.e. region is significant.


em=emmeans(bw_reg_fixed, ~"Reg")
pairs(em)

# contrast estimate     SE   df t.ratio p.value
# IM - LA  -0.05865 0.1160 2483  -0.505  0.9960
# IM - MG  -0.45431 0.1040 2483  -4.373  0.0002
# IM - NG  -0.05160 0.0944 2483  -0.546  0.9942
# IM - SG  -0.72288 0.1090 2483  -6.641  <.0001
# IM - SI   0.05518 0.1000 2483   0.550  0.9940
# LA - MG  -0.39567 0.1080 2483  -3.652  0.0036
# LA - NG   0.00705 0.1030 2483   0.068  1.0000
# LA - SG  -0.66424 0.1140 2483  -5.817  <.0001
# LA - SI   0.11382 0.1130 2483   1.012  0.9142
# MG - NG   0.40272 0.0879 2483   4.579  0.0001
# MG - SG  -0.26857 0.0924 2483  -2.907  0.0428
# MG - SI   0.50949 0.1010 2483   5.027  <.0001
# NG - SG  -0.67129 0.0980 2483  -6.847  <.0001
# NG - SI   0.10677 0.0958 2483   1.114  0.8757
# SG - SI   0.77806 0.1050 2483   7.404  <.0001

reg_pred_f=ggpredict(bw_reg_fixed, condition=c(MotherStatus='Milk'), terms = c("Reg","Sex[1]","AgeHrs[0]"))%>%
  mutate(x = fct_relevel(x, "SI", "IM", "LA", "NG","MG",  "SG"))

reg_pred_m=ggpredict(bw_reg_fixed, condition=c(MotherStatus='Milk'), terms = c("Reg","Sex[2]","AgeHrs[0]"))%>%
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

## plotting using ggregplot
inla_bw_plot=ggField(IM_spde, Mesh)+
  labs(fill = "Capture weight (kg) \n(as deviation from \nmean in SD units)")+
  theme_bw()+
  scale_fill_discrete_sequential(palette = "Burg", rev=FALSE)+
  theme(text = element_text(size = 18),
        legend.title=element_text(size=rel(0.7))) +
  annotate("segment", x = 1372, xend = 1382, y = 8000, yend = 8000, colour = "black", linewidth = 1) +
  annotate("text" ,x = 1377, y = 8001.5, label = "1km")+
  geom_point(data = bw_df%>%
               filter(!E<1355)%>%
               filter(!E>1385)%>%
               filter(!N<7997.5), aes(x = E, y = N), alpha=0.15) # Specify data
  


inla_bw_plot


birth_weight=bw_reg+inla_bw_plot+plot_annotation(tag_levels = 'A')
birth_weight

save(bw_reg,inla_bw_plot, file = "Deer_spatial_variation_ID/plots/bw_plots.RData")


# ggsave(birth_weight,
#        file = "Deer_spatial_variation_ID/plots/bw_spatial_var.png",
#        width = 11,
#        height = 10,
#        dpi=1000)
