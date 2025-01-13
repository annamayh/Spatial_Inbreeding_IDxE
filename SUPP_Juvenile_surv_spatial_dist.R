library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(patchwork)
library(emmeans)
library(INLA)
library(inlabru)
library(ggregplot)
library(colorspace)

setwd("/Volumes/Seagate_HD/")

surv_loc_df=read.table("Deer_spatial_variation_ID/survival_loc_2024_new_calves.txt", sep = ",", header = TRUE)%>%
  select(-MumFROH)%>%
  filter(!E>1385)%>%
  filter(!N<7997.5) 

head(surv_loc_df)

surv_loc_df$Code=as.factor(surv_loc_df$Code)
surv_loc_df$MumCode=as.factor(surv_loc_df$MumCode)
surv_loc_df$BirthYear=as.factor(surv_loc_df$BirthYear)
surv_loc_df$Sex=as.factor(surv_loc_df$Sex)
surv_loc_df$MotherStatus=as.factor(surv_loc_df$MotherStatus)
surv_loc_df$Reg=as.factor(surv_loc_df$Reg)
## base model of juvenile survival
suv_model_simple=glmmTMB(juvenile_survival~ 1+ Sex + MotherStatus + mum_age+mum_age_sq+Day_seq+FROH+BirthWt+
                           (1|BirthYear)+(1|MumCode), 
                         family=binomial, 
                         data=surv_loc_df, 
                         na.action = na.omit)


surv_reg_fixed=update(suv_model_simple, ~ . + Reg) ##just region as fixed effect
summary(surv_reg_fixed)

em=emmeans(surv_reg_fixed, ~"Reg")
pairs(em)


reg_pred_f=ggpredict(surv_reg_fixed, terms = c("Reg", "Sex[1]"))%>%
  mutate(x = fct_relevel(x, "SI", "IM", "LA", "NG","MG",  "SG"))
reg_pred_m=ggpredict(surv_reg_fixed, terms = c("Reg", "Sex[2]"))%>%
  mutate(x = fct_relevel(x, "SI", "IM", "LA", "NG","MG",  "SG"))
reg_pred=rbind(reg_pred_f,reg_pred_m)%>%
  arrange(x)%>%
  na.omit()

reg_pred

surv_reg=reg_pred%>%
  ggplot(aes(x=x, y=predicted, color=x, ymin=conf.low, ymax=conf.high, group=group))+
  geom_pointrange(linewidth=1, position = position_dodge(width=0.5))+
  theme_bw()+
  scale_color_manual(values = c("#f0a30a" ,"#a20025","#00aba9","chocolate1","#60a917", "#647687"))+
  labs(x="Spatial region", y="Predicted juvenile survival probability")+
  theme(text = element_text(size = 18),legend.position = "none")
surv_reg





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
              Day_seq=surv_loc_df$Day_seq, 
              BirthWt=surv_loc_df$BirthWt)
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


IM_spde  <- inla(y~ -1 + Intercept+Sex + MotherStatus + mum_age+mum_age_sq+Day_seq+FROH+BirthWt+
                   f(BirthYear, model = 'iid')+f(MumCode, model = 'iid')+ f(w, model=spde), 
                 family = "binomial",
                 data=inla.stack.data(stackfit), 
                 control.compute = list(dic=TRUE),
                 control.predictor = list(
                   A=inla.stack.A(stackfit))) 


summary(IM_spde)



inla_surv_plot=ggField(IM_spde, Mesh)+
  labs(fill = "Juvenile survival \n(as untransformed \ndevaition from mean)")+
  theme_bw()+
  scale_fill_discrete_sequential(palette = "Oranges", rev=FALSE)+
  theme(text = element_text(size = 18),
        legend.title=element_text(size=rel(0.7))) +
  annotate("segment", x = 1372, xend = 1382, y = 8000, yend = 8000, colour = "black", linewidth = 1) +
  annotate("text" ,x = 1377, y = 8001.5, label = "1km")+
  geom_point(data = surv_loc_df%>%
               filter(!E<1355)%>%
               filter(!E>1385)%>%
               filter(!N<7997.5), aes(x = E, y = N), alpha=0.15) # Specify data

inla_surv_plot



supp_surv=surv_reg+inla_surv_plot+plot_annotation(tag_levels = 'A')

ggsave(supp_surv, 
       file = "Deer_spatial_variation_ID/plots/SUPP_juve_surv_wBW_updated.png",
       width = 11,
       height = 6,
       dpi=1000
       )

##updated 27.11.24
