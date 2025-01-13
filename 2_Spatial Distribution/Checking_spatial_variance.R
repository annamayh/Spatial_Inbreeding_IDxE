library(tidyverse)

setwd("/Volumes/Seagate_HD/") #external hardrive

# mum status at conception
mum_stat=read.csv("Deer_spatial_variation_ID/Rum_data_files/sys_HindStatusAtConception.csv") %>%
  select(-CalfBirthYear)%>%rename(MumCode=Mum, Code=Calf)
#census records
census_full=read.csv("Deer_spatial_variation_ID/Rum_data_files/tblCensus_Nov2024.csv")%>%
  select(Date,Code,Northing,Easting) 
# birth year of mother to get mothers age at conception
mum_birthyr=read.csv("Deer_spatial_variation_ID/Rum_data_files/tbllife.csv")%>%
  select(Code,BirthYear)%>%rename(mum_birthyear=BirthYear,MumCode=Code) 
# birth weight (estimated and capture weight)
birth_wt<-read.csv("Deer_spatial_variation_ID/Rum_data_files/sys_BirthWt.csv")
# death and life info
life=read.csv("Deer_spatial_variation_ID/Rum_data_files/tbllife.csv")%>%
  select(Code, BirthDay, BirthMonth, BirthYear, DeathDay, DeathMonth, DeathYear, Sex, MumCode, DeathType)
# last recorded sighting
last_cen_sight=read.csv("Deer_spatial_variation_ID/Rum_data_files/sys_LastCensusSighting.csv")


####################################################################################
## Getting the average location of the mother in calves firs year ##
####################################################################################

census_yr=census_full%>%
  separate_wider_delim(Date,names=c("Year","Month","Day"),delim="-", cols_remove = F)%>%
  mutate(
    Month=as.numeric(Month),
    Year=as.numeric(Year),
    deer_year = ifelse(Month < 5, Year - 1, Year) ## deer year runs from from 1st May to 30th April, 
    # so if month is <5 then deer year is assigned to previous yr
  )

ave_N_E_deer_year=census_yr%>%
  filter(Code!='')%>%
  group_by(Code, deer_year)%>%
  summarise(
    N_mean = mean(Northing, na.rm = TRUE), # getting the mean northing and easting for each deer year
    E_mean = mean(Easting, na.rm = TRUE), 
    N_SD = sd(Northing), 
    E_SD = sd(Easting)
  ) %>%
  ungroup() %>%
  na.omit()

#average location of the mum during the birth year of the calf
mum_loc_calf_yr=ave_N_E_deer_year%>%
  rename(MumCode=Code, BirthYear=deer_year)%>% 
  inner_join(life, by = c('BirthYear', 'MumCode'))%>%
  select(Code, MumCode, BirthYear, N_mean, E_mean, N_SD, E_SD)%>%
  ## assigning them a region based on ave N and E 
  mutate(Reg = case_when(
    N_mean < 8019 ~ "SG", # South Glen
    E_mean < 1361 ~ "LA", # Laundry Greens
    E_mean < 1366 & N_mean > 8033 ~ "NG", # North Glen
    E_mean < 1366 & N_mean <= 8033 ~ "MG", # Mid Glen
    E_mean < 1373 ~ "IM", # IM
    TRUE ~ "SI" # SI
  ))


## all used data
dat=read.table("Deer_spatial_variation_ID/df_loc_2024_new_calves_all.txt", sep = ",", header = TRUE)


dat_w_var=left_join(dat, mum_loc_calf_yr, by = join_by(Code, BirthYear, MumCode, Reg))
  


mean(dat_w_var$N_SD,na.rm = T)
median(dat_w_var$N_SD, na.rm = T)

max(dat_w_var$N_SD, na.rm = T)
min(dat_w_var$N_SD, na.rm = T)

mean(dat_w_var$E_SD, na.rm = T)
median(dat_w_var$E_SD, na.rm = T)
max(dat_w_var$E_SD, na.rm = T)
min(dat_w_var$E_SD, na.rm = T)

nrow(filter(dat_w_var$N_SD))
