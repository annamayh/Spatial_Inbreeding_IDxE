library(tidyverse)
library(stringi)

setwd("/Volumes/Seagate_HD/")

all_anon=read.table("Deer_spatial_variation_ID/df_loc_2024_new_calves_all.txt", sep = ",", header = TRUE)

ids=all_anon$Code # all individual ids
mum_ids=unique(all_anon$MumCode) # all mum codes


all_names=union(ids, mum_ids) # get a list of both (removing duplicates)

# randomly generate new names for all of these ids
new_names <- paste0(
  stri_rand_strings(length(all_names), 4, pattern = "[A-Z]"), 
  stri_rand_strings(length(all_names), 1, pattern = "[0-9]")
)

# link new names to old names
id_mapping <- setNames(new_names, all_names)

# Replace Codes with totally randomised codes 
all_anon_fin <- all_anon %>%
  mutate(
    Code = id_mapping[Code],
    MumCode = id_mapping[MumCode]
  )



write.table(all_anon_fin,
            file = "Deer_spatial_variation_ID/Full_dataset_IDxE_anon.txt",
            row.names = F, quote = F, sep = ",",na = "NA") #saving tables as txt file 

