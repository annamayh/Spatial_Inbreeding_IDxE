library(tidyverse)

setwd("/Volumes/Seagate Por") #external hardrive

# mum status at conception
mum_stat=read.csv("Deer_spatial_variation_ID/Rum_data_files/sys_HindStatusAtConception.csv") %>%
  select(-CalfBirthYear)%>%rename(MumCode=Mum, Code=Calf)
#census records
census_full=read.csv("Deer_spatial_variation_ID/Rum_data_files/tblCensus.csv")%>%
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
  separate_wider_delim(Date,names=c("Day","Month","Year"),delim="/", cols_remove = F)%>%
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
    E_mean = mean(Easting, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  na.omit()

#average location of the mum during the birth year of the calf
mum_loc_calf_yr=ave_N_E_deer_year%>%
  rename(MumCode=Code, BirthYear=deer_year)%>% 
  inner_join(life, by = c('BirthYear', 'MumCode'))%>%
  select(Code, MumCode, BirthYear, N_mean, E_mean)%>%
  ## assigning them a region based on ave N and E 
  mutate(Reg = case_when(
    N_mean < 8019 ~ "SG", # South Glen
    E_mean < 1361 ~ "LA", # Laundry Greens
    E_mean < 1366 & N_mean > 8033 ~ "NG", # North Glen
    E_mean < 1366 & N_mean <= 8033 ~ "MG", # Mid Glen
    E_mean < 1373 ~ "IM", # IM
    TRUE ~ "SI" # SI
  ))

## check regions make sense
ggplot(mum_loc_calf_yr, aes(y=N_mean, x=E_mean,color=Reg))+
  geom_point()+
  theme_classic()


###########################################
#### Calculate juvenile survival #####
#########################################

census_sep=last_cen_sight%>%
  mutate(LastSeen_date=as.Date(LastSeen, format="%d.%b.%y"))%>%
  separate(LastSeen_date,into=c("YearLast","MonthLast","DayLast"),sep="-")%>%
  mutate(YearLast=as.numeric(YearLast))%>%
  select(Code, YearLast)%>%
  left_join(life)%>%
  mutate(yrs_since_last_Seen=(2024-YearLast))

## for all ids with missing death years fill in estimated death year when havent been seen for at least 3 yrs 
year_plus_census_filt=census_sep%>%
  filter(is.na(DeathYear))%>%
  mutate(DeathYear_est=YearLast+3)%>%
  filter(yrs_since_last_Seen>1)%>% #has to be more than 3 yrs scince its been seen to have an estimated death
  filter(DeathYear_est<2024)%>% 
  select(Code,DeathYear_est,YearLast)

##adding in an estimated death yr and year last seen for ids w/ no death year (not including ids born in the last 3 yrs )
year_plus_est_death=life%>%left_join(year_plus_census_filt) %>%
  mutate(
    DeathYear=as.numeric(DeathYear),
    DeathYear_est=as.numeric(DeathYear_est),
    BirthYear=as.numeric(BirthYear),
    YearLast=as.numeric(YearLast)
  )

year_juve=year_plus_est_death%>%
  mutate(juvenile_survival = case_when(
    ((is.na(DeathYear)) & (BirthYear==YearLast))~"0",
    
    ## when there is no recorded death year use the estimated death year (which is when id has not been seen for 3 years...)
    # not sure any juveniles will actually be recorded as 0 in this case as theyre killed off after 3 year anyway 
    (is.na(DeathYear) & (DeathYear_est<(BirthYear+2)))~"0",
    #if died in same year born then juvenile surv=0
    (DeathYear==BirthYear) ~ "0", 
    #if died the year after it was born juvenle surv=0
    (DeathYear==(BirthYear+1)) ~ "0", 
    #if died before May in the second year of life then did not survival to juvenile (cut off is May)
    (DeathYear==(BirthYear+2))&(DeathMonth<5) ~ "0", 
    #any ids with no recorded death month and a death year 2 yrs after birth is given NA as we dont know if id died before or after 2nd birhday (only 5 ids)
    (DeathYear==(BirthYear+2))&(is.na(DeathMonth)) ~ "NA", 
    
  ) )%>%
  # all ids that didnt get a 0 for juvenile survival gets a 1 for successfully surviving 
  mutate(juvenile_survival = replace_na(juvenile_survival,"1"))


## filtering out shot and accidental deaths (and those that are still alive)
juvenile_surv=year_juve%>%filter(DeathType!= "S"| is.na(DeathType))%>% #filter out ids that were shot
  filter(DeathType!= "A" | is.na(DeathType)) %>% ##filter out ids that were an accidental death
  filter(DeathType!= "D" | is.na(DeathType))%>% ## keeping NA death type
  select(Code, BirthYear, Sex, MumCode, juvenile_survival)

table(juvenile_surv$juvenile_survival)

# want to keep ids that were shot but also made it to adulthood 
juvenile_surv_shot=year_juve%>%filter(DeathType== "S"&juvenile_survival=="1")%>%
  select(Code, BirthYear, Sex, MumCode, juvenile_survival)

table(juvenile_surv_shot$juvenile_survival)

## adding ids that survived to age 2 then were shot to the df
juvenile_surv=juvenile_surv%>%rbind(juvenile_surv_shot)
table(juvenile_surv$juvenile_survival)


##############################################
####    DOB continuous   ####
###############################################
# first birth date is 30th/4
# Last is 10th Oct
day_seq=read.csv("Deer_spatial_variation_ID/Rum_data_files/date_seq.csv", header = TRUE)%>%
  mutate(year_date=as.Date(Date, format="%d.%b"))%>%
  separate(year_date,into=c(NA,"BirthMonth","BirthDay"),sep="-")%>%
  mutate(
    BirthDay=as.numeric(BirthDay),
    BirthMonth=as.numeric(BirthMonth))

DOBs=year_juve%>%select(Code, BirthDay, BirthMonth)%>%
  left_join(day_seq, by = c('BirthDay', 'BirthMonth'))%>% ## sequences of day starting at 30th April 
  select(Code, Day_seq)


##############################################
## FROH values ##
##############################################
FROH_full<-read.table("PhD/2023_ROH_search/2021_sleuthed_052023.hom.indiv", header=T, stringsAsFactors = F)%>%
  select(IID,KB) %>% rename(Code=IID)%>%mutate(FROH=KB/2591865)%>%
  filter(nchar(Code)==5)%>% #removing IDs with non-sensical ID codes
  select(-KB)

FROH_mum<-read.table("PhD/2023_ROH_search/2021_sleuthed_052023.hom.indiv", header=T, stringsAsFactors = F)%>%
 select(IID,KB) %>% rename(MumCode=IID)%>%mutate(MumFROH=KB/2591865)%>%
  filter(nchar(MumCode)==5)%>% #removing IDs with non-sensical ID codes
  select(-KB)

mum_loc_calf_yr$BirthYear=as.numeric(as.character(mum_loc_calf_yr$BirthYear))


juvenile_loc_surv_bw=juvenile_surv%>%
  left_join(mum_stat)%>%
  left_join(birth_wt)%>%
  left_join(mum_birthyr)%>%
  mutate(mum_age=BirthYear-mum_birthyear)%>%
  mutate(mum_age_sq=mum_age^2)%>%
  left_join(mum_loc_calf_yr)%>%
  left_join(FROH_full)%>%
  left_join(FROH_mum)%>%
  rename(N=N_mean, E=E_mean)%>%
  left_join(DOBs)
  
#selecting columns for juvenile survival 
surv_loc_df= juvenile_loc_surv_bw%>%
  select(-mum_birthyear, -CaptureWt, -AgeHrs) %>%
  drop_na(FROH, juvenile_survival)

head(surv_loc_df)
nrow(surv_loc_df)

## selecting columnc for birth weight (uses capture weight and age in hrs )
bw_loc_df= juvenile_loc_surv_bw%>%
  select(-mum_birthyear, -BirthWt,-juvenile_survival) %>%
  drop_na(FROH,CaptureWt)

head(bw_loc_df)
nrow(bw_loc_df)



write.table(surv_loc_df,
            file = "Deer_spatial_variation_ID/survival_loc_2024.txt",
            row.names = F, quote = F, sep = ",",na = "NA") #saving tables as txt file 


write.table(bw_loc_df,
            file = "Deer_spatial_variation_ID/birthweight_loc_2024.txt",
            row.names = F, quote = F, sep = ",",na = "NA") #saving tables as txt file 

