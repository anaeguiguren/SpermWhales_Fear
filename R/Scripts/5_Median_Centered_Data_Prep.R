# Median Centered Data prep------
rm(list = ls())
source("R/Scripts/2_Prep_Data_For_Analysis.R")

library(factoextra)
library(scales)
library(dplyr)
library(patchwork)
library(ggforce)
library(lubridate)

#1. GET DAILY MEDIANS----

day_medians <- dat %>%
  group_by(date) %>%
  summarize(
    NumW = median(NumW, na.rm = T),
    sd_head = median(sd_head, na.rm = T),
    Rank = median(Rank, na.rm = T),
    File = median(File, na.rm = T),
    PVA = median(PVA, na.rm = T),
    Mean_NN = median(Mean_NN, na.rm = T), 
    mean_body_angle = median(mean_body_angle, na.rm = T),
    Maj_axis = median(Maj_axis, na.rm = T),
    soc_med = median(soc, na.rm = T),
    n = n()
  )



# 2. CALCULATE MEDIAN CENTERED VALUES for full dataset-----

#create identifying index
idx <- match(dat$date, day_medians$date)

#name variables to modify
vars <- c("NumW", "sd_head", "Rank", "File", "PVA", "Mean_NN", 
          "mean_body_angle", "Maj_axis")

#create new data frame for median centered data

dat_med_cent <- dat
dat_med_cent[vars] <- dat_med_cent[vars] - day_medians[idx, vars]

#make useful behav state
dat_med_cent <- 
  dat_med_cent %>%
  mutate(behav_state_int = ifelse(behav_state == "for/brea/tra", 
                              "not_interacting", 
                              ifelse(behav_state == "socializing", 
                                     "not_interacting", 
                                     as.character(behav_state))))



# 3. TIDY AND SAVE-----
# keep relevant columns

dat_med_cent <- dat_med_cent %>%
  select(date,flight_session_ID, frame_local_time, frame_ID, group_ID,
         event_state,
         behav_state_int,
         intensity_sp, 
         NumW, sd_head, Rank, File, PVA, Mean_NN, 
         mean_body_angle, Maj_axis, soc_raw = soc
         )


# save 
saveRDS(dat_med_cent, 
        "Data/Processed_Data/median_centered_full_data_form.rds")


#save d for processing
saveRDS(d, 
        "Data/Processed_Data/event_descriptions.rds")
