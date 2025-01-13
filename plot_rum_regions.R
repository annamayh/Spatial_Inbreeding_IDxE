library(tidyverse)

setwd("/Volumes/Seagate_HD/")


rum_outline=read.csv("PhD/Chapter_5_spatial_ID_x_E/Spatial_var_inbreeding/INLA/RumBoundary.csv")%>%
  rename(E = Easting, N = Northing) 

surv_loc_df=read.table("Deer_spatial_variation_ID/df_loc_2024_new_calves_all.txt", sep = ",", header = TRUE)%>%
  filter(!E<1355)%>%
  filter(!E>1385)%>%
  filter(!N<7997.5) %>%#remove 38 records outside the syudy areas
  drop_na(FROH)

surv_loc_df%>%
  group_by(Reg)%>%
  summarise(n=n())

cols = c("SI"="#f0a30a" ,"IM"="#a20025","LA"="#00aba9","NG"="chocolate1","MG"="#60a917", "SG"="#647687")
reg_labs= c('SI'="Shamhnan Insir (SI)\n n = 558", "IM"="Intermediate area (IM)\n n = 441", "LA"="Laundry greens (LA)\n n = 439", "NG"="North glen (NG)\n n = 615", "MG"="Mid glen (MG)\n n = 389", "SG"="South glen (SG)\n n = 348")

rum_map=ggplot() +
  geom_polygon(data = rum_outline, aes(x = E, y = N), alpha = 0.25, colour='black') +  # Specify data and group
  geom_point(data = surv_loc_df, aes(x = E, y = N, group = Reg, colour = Reg), alpha=0.6) + # Specify data
  scale_color_manual(values=cols, labels=reg_labs, breaks=c("SI","IM","LA", "NG", "MG", "SG")) +
  theme_classic() +
  theme(text = element_text(size = 18)) +
  labs (y="Northing", x='Easting', colour = 'Spatial region')+
  annotate("segment", x = 1372, xend = 1382, y = 8000, yend = 8000, colour = "black", linewidth = 1) +
  annotate("text" ,x = 1377, y = 8001.5, label = "1km")

rum_map

ggsave(rum_map,
       file = "Deer_spatial_variation_ID/REVISED/Figs/Fig1.jpeg",
       width = 6.5,
       height = 5, 
       bg = "white"
)
# updated 25.11.23