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

#mupdated 25.11.23 ## 
# to transform FROH using sqrt # 

setwd("/Volumes/Seagate_HD/")

surv_loc_df=read.table("Deer_spatial_variation_ID/df_loc_2024_new_calves_all.txt", sep = ",", header = TRUE)%>%
  filter(!E>1385)%>%
  filter(!N<7997.5)%>%#remove 38 records outside the syudy areas
  drop_na(FROH)%>%
  mutate(FROH_trans=sqrt(FROH))

mean(surv_loc_df$FROH)
sd(surv_loc_df$FROH)


surv_loc_df%>%
  group_by(Reg)%>%
  summarise(mean=mean(FROH))


##############################################################################
#### testing difference region as categorical ##############################
###########################################################################
surv_loc_df$BirthYear=as.numeric(surv_loc_df$BirthYear)


FROH_reg=surv_loc_df%>%
  select(Code, BirthYear, Sex, MumCode, FROH_trans, Reg, N, E)%>%
  filter(Sex!=3)%>%
  na.omit()%>%
  mutate(year_cont=BirthYear-min(BirthYear))#



FROH_reg$Code=as.factor(FROH_reg$Code)
FROH_reg$MumCode=as.factor(FROH_reg$MumCode)
FROH_reg$BirthYear=as.factor(FROH_reg$BirthYear)
FROH_reg$Sex=as.factor(FROH_reg$Sex)
FROH_reg$Reg=as.factor(FROH_reg$Reg)


cohen.d(FROH_reg$FROH[FROH_reg$Reg == "SI"], FROH_reg$FROH[FROH_reg$Reg == "SG"])
#cohen.d(FROH_reg$FROH[FROH_reg$Reg == "SI"], FROH_reg$FROH[FROH_reg$Reg == "LA"])
#cohen.d(FROH_reg$FROH[FROH_reg$Reg == "SI"], FROH_reg$FROH[FROH_reg$Reg == "NG"])
cohen.d(FROH_reg$FROH[FROH_reg$Reg == "SI"], FROH_reg$FROH[FROH_reg$Reg == "NG"])


cohen.d(FROH_reg$FROH[FROH_reg$Reg == "IM"], FROH_reg$FROH[FROH_reg$Reg == "SG"])
#cohen.d(FROH_reg$FROH[FROH_reg$Reg == "IM"], FROH_reg$FROH[FROH_reg$Reg == "LA"])
cohen.d(FROH_reg$FROH[FROH_reg$Reg == "IM"], FROH_reg$FROH[FROH_reg$Reg == "NG"])
#cohen.d(FROH_reg$FROH[FROH_reg$Reg == "IM"], FROH_reg$FROH[FROH_reg$Reg == "MG"])



FROH_model_base=glmmTMB(FROH_trans~ year_cont+
                            (1|MumCode), 
                          family=gaussian(), 
                          data=FROH_reg, 
                          na.action = na.omit,
)

FROH_model_simple=update(FROH_model_base, ~ . + Reg) ##just region as fixed effect

anova(FROH_model_base,FROH_model_simple) ## region fixed is a better fit i.e. region is significant.


summary(FROH_model_simple)

em=emmeans(FROH_model_simple, ~"Reg")
pairs(em)

#updated 25.11.23
# 
# contrast estimate      SE   df t.ratio p.value
# IM - LA   0.00179 0.00379 2779   0.471  0.9971
# IM - MG   0.00760 0.00378 2779   2.011  0.3364
# IM - NG   0.01056 0.00340 2779   3.104  0.0237 <<
# IM - SG   0.01538 0.00387 2779   3.975  0.0010 <<
# IM - SI  -0.00203 0.00349 2779  -0.582  0.9922
# LA - MG   0.00581 0.00380 2779   1.528  0.6463
# LA - NG   0.00878 0.00348 2779   2.524  0.1174
# LA - SG   0.01359 0.00391 2779   3.479  0.0068 <<
# LA - SI  -0.00381 0.00358 2779  -1.067  0.8945
# MG - NG   0.00296 0.00344 2779   0.862  0.9554
# MG - SG   0.00778 0.00383 2779   2.028  0.3265
# MG - SI  -0.00963 0.00355 2779  -2.709  0.0738 <<
# NG - SG   0.00481 0.00358 2779   1.346  0.7590
# NG - SI  -0.01259 0.00319 2779  -3.949  0.0011 <<
# SG - SI  -0.01741 0.00365 2779  -4.763  <.0001 <<

#P value adjustment: tukey method for comparing a family of 6 estimates 

FROH_reg_pred=ggpredict(FROH_model_simple, terms = c("Reg"))

reg=FROH_reg_pred%>%
  mutate(x = fct_relevel(x, "SI", "IM", "LA", "NG","MG","SG" ))%>%
  mutate(predicted_rescaled = predicted^2, 
         conf.low_rescaled = conf.low^2, 
         conf.high_resscaled = conf.high^2
         )%>%
  na.omit()%>%
  ggplot( aes(x=x, y=predicted_rescaled, color=x, ymin=conf.low_rescaled, ymax=conf.high_resscaled))+
  geom_pointrange(linewidth=1)+
  theme_bw()+
  scale_color_manual(values = c("#f0a30a" ,"#a20025","#00aba9","chocolate1","#60a917", "#647687"))+
  labs(x="Spatial region", y=expression(paste("Predicted F"[ROH])))+
  theme(text = element_text(size = 18),legend.position = "none")
reg





########################################################################
############# region as matrix INLA ##################################
########################################################################
## simple models first
IM1  <- inla(FROH_trans~year_cont, 
             family = "gaussian",
             data = FROH_reg,
             control.compute = list(dic=TRUE)) 

IM2  <- inla(FROH_trans~year_cont+ f(MumCode, model = 'iid'), 
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
  data=list(y=FROH_reg$FROH_trans), 
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
  labs(fill = "FROH\n(deviation \nfrom mean \nin SD units)")+
  theme_bw()+
  scale_fill_continuous_sequential(palette = "BluYl")+
  theme(text = element_text(size = 18),
        legend.title=element_text(size=rel(0.8))) +
  annotate("segment", x = 1372, xend = 1382, y = 8000, yend = 8000, colour = "black", linewidth = 1) +
  annotate("text" ,x = 1377, y = 8001.5, label = "1km")+
  geom_point(data = FROH_reg%>%
               filter(!E<1355)%>%
               filter(!E>1385)%>%
               filter(!N<7997.5), aes(x = E, y = N), alpha=0.1) # Specify data
# havent transformed these because they are deviations from mean anyway



froh_spat=reg+inla_froh_gg+plot_annotation(tag_levels = 'A')
froh_spat                    

ggsave(froh_spat,
       file = "Deer_spatial_variation_ID/REVISED/Figs/Fig2.jpeg",
       width = 9,
       height = 5, 
       bg = "white"
) #updated 25.11.23

