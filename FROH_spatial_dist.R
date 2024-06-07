library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(patchwork)
library(emmeans)
library(INLA)
library(inlabru)
library(ggregplot)
library(colorspace)
library(RColorBrewer)


setwd("/Volumes/Seagate Por")

surv_loc_df=read.table("PhD/Chapter_5_spatial_ID_x_E/Spatial_var_inbreeding/survival_loc.txt", sep = ",", header = TRUE)%>%
  filter(!E>1385)%>%
  filter(!N<7997.5)#remove 38 records outside the syudy areas


## FROH per region violin plots 

cols <- c("SI"="#f0a30a" ,"IM"="#a20025","LA"="#00aba9","NG"="chocolate1","MG"="#60a917", "SG"="#647687")

ggplot(surv_loc_df, aes(x=Reg, y=FROH, colour=Reg))+
  scale_x_discrete(limits=c("SI", "IM", "LA", "NG","MG",  "SG"))+
  geom_violin()+
  coord_flip()+
  geom_jitter(position=position_jitter(0.2), alpha=0.5)+
  scale_color_manual(values=cols)+
  theme_bw()+
  geom_hline(yintercept = 0.1, linetype=2)+
  geom_hline(yintercept = 0.2, linetype=2)+
  geom_hline(yintercept = 0.3, linetype=2)+
  theme(text = element_text(size = 18), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")+
  labs(y=expression(F["ROH"]))
  

##############################################################################
#### testing difference region as categorical ##############################
###########################################################################
surv_loc_df$BirthYear=as.numeric(surv_loc_df$BirthYear)


FROH_reg=surv_loc_df%>%
  select(Code, BirthYear, Sex, MumCode, FROH, Reg, N, E)%>%
  filter(Sex!=3)%>%
  na.omit()%>%
  mutate(year_cont=BirthYear-min(BirthYear))#



FROH_reg$Code=as.factor(FROH_reg$Code)
FROH_reg$MumCode=as.factor(FROH_reg$MumCode)
FROH_reg$BirthYear=as.factor(FROH_reg$BirthYear)
FROH_reg$Sex=as.factor(FROH_reg$Sex)
FROH_reg$Reg=as.factor(FROH_reg$Reg)



FROH_model_simple=glmmTMB(FROH~ year_cont+ Reg+
                            (1|MumCode), 
                          family=gaussian(), 
                          data=FROH_reg, 
                          na.action = na.omit,
)



summary(FROH_model_simple)


FROH_reg_pred=ggpredict(FROH_model_simple, terms = c("Reg"))

reg=FROH_reg_pred%>%
  mutate(x = fct_relevel(x, "SI", "IM", "LA", "NG","MG",  "SG"))%>%
  ggplot( aes(x=x, y=predicted, color=x, ymin=conf.low, ymax=conf.high))+
  geom_pointrange(linewidth=1)+
  theme_bw()+
  scale_color_manual(values = c("#f0a30a" ,"#a20025","#00aba9","chocolate1","#60a917", "#647687"))+
  labs(x="Spatial region", y=expression(paste("F"[ROH])))+
  theme(text = element_text(size = 18),legend.position = "none")
reg

 

########################################################################
############# region as matrix INLA ##################################
########################################################################
## simple models first
IM1  <- inla(FROH~year_cont, 
             family = "gaussian",
             data = FROH_reg,
             control.compute = list(dic=TRUE)) 

IM2  <- inla(FROH~year_cont+ f(MumCode, model = 'iid'), 
             family = "gaussian",
             data = FROH_reg,
             control.compute = list(dic=TRUE)) 

## read in rum
rum_outline=read.csv("PhD/Chapter_5_spatial_ID_x_E/Spatial_var_inbreeding/INLA/RumBoundary.csv")%>%
  rename(E = Easting, N = Northing) 


N=nrow(rum_outline)
rum_line_rev=rum_outline[N:1, c("E","N")]

Mesh=inla.mesh.2d(loc.domain= rum_outline, 
                  max.edge=2, #probs use 1 for actual model
                  boundary=
                    inla.mesh.segment(rum_line_rev))

plot(Mesh, asp=1)


#################################################################################################
## setting up INLA model ########################################################################
######################################################################################################

#set weighting using A matrix
Locations=cbind(FROH_reg$E, FROH_reg$N)#locations of ids

A=inla.spde.make.A(Mesh, loc=Locations)

#define SPDE 
spde=inla.spde2.matern(Mesh, alpha = 2) # would need to adjust for time series data

#define spatial field
w.index=inla.spde.make.index(name = 'w', n.spde=spde$n.spde, n.group = 1, n.repl = 1)

#make model matrix
N <- nrow(FROH_reg)
X0=data.frame(Intercept = rep(1, N),
              year_cont = FROH_reg$year_cont)

x=as.data.frame(X0)


## now fitting spde and birth year as random effects 
# have to re-do the stack
stackfit2=inla.stack(
  data=list(y=FROH_reg$FROH), 
  A = list(1, 1, A), 
  effects=list(
    X=x,
    MumCode = FROH_reg$MumCode, # insert vectors of any random effects
    w=w.index)
)


Fspat2=as.formula(paste0("y ~ -1 + Intercept + year_cont + f(MumCode, model = 'iid') + f(w, model=spde) "))

#run model
IM_sp2=inla(y ~ -1 + Intercept + year_cont + 
              f(MumCode, model = 'iid') + f(w, model=spde),
            family = "gaussian", 
            data=inla.stack.data(stackfit2), 
            control.compute = list(dic=TRUE),
            control.predictor = list(
              A=inla.stack.A(stackfit2))#, verbose = TRUE
            
)

summary(IM_sp2)




## check whether DIC decreases with SPDE 
SpatialList <- list(IM1, IM2, IM_sp2)
sapply(SpatialList, function(f) f$dic$dic)
INLADICFig(SpatialList, ModelNames = c("IM1","IM2" ,"SPDE_1"))

library(viridisLite)
coul <- viridis(100)


inla_froh_gg=ggField(IM_sp2, Mesh, Fill="Continuous")+
  labs(fill = "FROH\n(deviation \nfrom mean)")+
  theme_bw()+
  scale_fill_continuous_sequential(palette = "BluYl")+
  theme(text = element_text(size = 18),
        legend.title=element_text(size=rel(0.8)))


reg+inla_froh_gg+plot_annotation(tag_levels = 'A')
                     
