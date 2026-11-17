# GLMM on PCA-R outcomes
# Takes PCA-R outcomes, filters only frames observed within 2 hrs of event
# runs GLMM selection and plots effect sizes

#0. LOAD LIBRARIES-----
rm(list = ls())

#list of files
packages <- c("MuMIn", "ggplot2", "dplyr", "tidyr", 
              "lubridate", 
              "nlme","emmeans", "here")

#check if installed and install if needed

installed <- rownames(installed.packages())

for (p in packages) {
  if (!p %in% installed) {
    install.packages(p)
  }
  library(p, character.only = TRUE)
}

#load

library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)
library(here)

#1. LOAD RAW DATA-----
source("R/Scripts/2_Prep_Data_For_Analysis.R")

d <- readRDS("Data/Processed_Data/event_descriptions.rds")

#2. FILTER FOR EVENT BEFORE/AFTER-----
# keep only observations collected during interactions or
# up to 2 hours before/after, excluding anything within 10 mins
# of the start/end of an interaction

# format date time
dat <- dat %>%
  mutate(frame_local_time_f = as.POSIXct(frame_local_time, format = "%Y-%m-%d %H:%M:%S")) %>%
  arrange(frame_local_time_f)


#~2a. Get event information-----
#load data
e_raw <- read.csv("Data/Processed_Data/Boris_Event_Times.csv", header = T)

#get start and end of each interaction (focal follow based)
e <- e_raw %>%
  group_by(event) %>%
  summarize(start_date_time = first(start_date_time), 
            end_date_time = last(end_date_time), 
            duration_focal = sum(duration)) %>%
  mutate(start_date_time = as.POSIXct(start_date_time, format = "%Y-%m-%d %H:%M:%S"), 
         end_date_time = as.POSIXct(end_date_time, format = "%Y-%m-%d %H:%M:%S") )

# fix event names for matching format

e <- e %>%
  mutate(event.code = tolower(event))

# join with event descriptions 

e <- left_join(e, d, by = 'event.code')

#~2b. create selection window-----

# make a selection window:
event_windows <- e %>% 
  mutate(
    # look at frames collected <2 hrs from events 
    start_include = start_date_time - hours(2), 
    end_include = end_date_time + hours(2), 
    
    # exclude frames collected within 10 mins of start/end events
    start_exclude = start_date_time - minutes(10), 
    end_exclude = end_date_time + minutes(10)
    
  )


#3. ASSIGN WINDOW INFO TO DATA----

#~3a. join dat with event windows to get event times----
#for each frame.
dat_with_events <- dat %>%
  
  # Cross join
  left_join(event_windows, by = character()) %>%
  
  # Calculate timing metrics
  mutate(
    time_to_start = as.numeric(difftime(frame_local_time_f,
                                        start_date_time,
                                        units = "secs")),
    
    time_to_end = as.numeric(difftime(frame_local_time_f,
                                      end_date_time,
                                      units = "secs")),
    
    abs_time_to_event = pmin(abs(time_to_start),
                             abs(time_to_end),
                             na.rm = TRUE),
    
    # Same calendar day?
    same_day = as.Date(frame_local_time_f) ==
      as.Date(start_date_time),
    
    # Within extended event window?
    in_window = frame_local_time_f >= start_include &
      frame_local_time_f <= end_include
  ) %>%
  
  # KEEP ONLY reasonable candidate events
  filter(same_day | in_window) %>%
  
  # Choose closest candidate event
  group_by(frame_local_time_f) %>%
  slice_min(abs_time_to_event,
            n = 1,
            with_ties = FALSE) %>%
  ungroup()





#~3b. Assign interaction_window label----
dat_with_events <- dat_with_events %>%
  mutate(
    interaction_window = case_when(
      !is.na(start_date_time) & frame_local_time_f >= start_date_time & frame_local_time_f <= end_date_time ~ "during",
      !is.na(start_exclude) & frame_local_time_f >= start_exclude & frame_local_time_f <= end_exclude &
        (frame_local_time_f < start_date_time | frame_local_time_f > end_date_time) ~ "within 10 mins",
      !is.na(start_include) & frame_local_time_f >= start_include & frame_local_time_f <= end_include &
        (frame_local_time_f < start_exclude | frame_local_time_f > end_exclude) ~ "within 2hrs",
      as.Date(frame_local_time_f) == as.Date(start_date_time) ~ "more than 2 hours",
      TRUE ~ NA_character_
    )
  )

#~3c. assign a phase (before/after)----

# label as before/during/after
dat_with_events <- dat_with_events %>%
  mutate(
    phase = case_when(
      frame_local_time_f >= start_date_time & frame_local_time_f <= end_date_time ~ "during",
      frame_local_time_f <  start_date_time ~ "before",
      frame_local_time_f >  end_date_time   ~ "after"
    ), 
    event.code = event.code.y
  )


#~3d. keep only events with at least two phases including during----

# get summary of events by phase
events <- dat_with_events %>%
  group_by(event.code, phase) %>%
  summarise(n = n(), 
            date = first(date), 
            .groups = "drop") %>%
  group_by(event.code)%>%
  filter(sum(n) > 10 ) %>%
  ungroup() %>%
  arrange(as.Date(date))
  
  



# identify useful events
evs <- events %>%
  group_by(event.code) %>% 
  filter(any(phase == "during") & n_distinct(phase)>1) %>%
  select(event.code)%>%
  distinct(event.code)


dat_with_events <- dat_with_events %>%
  filter(event.code %in% evs$event.code)

#remove frames within 10 minutes
dat_with_events <- dat_with_events %>%
  filter(interaction_window != "within 10 mins" |
           is.na(interaction_window))

#~4e. cleanup variable names ------
dat_with_events <- dat_with_events %>%
  mutate(intensity_sp = intensity_sp.y, 
         intensity = intensity_sp.y,
         treat = interaction_window)





#4. MAKE SUMMARY TABLE-----
#re-sort factor names 
dat_with_events <- dat_with_events %>%
  mutate(phase = factor(phase, levels = c("before", "during", "after")), 
         treat = factor(treat, levels = c("more than 2 hours",
                                          "within 2hrs",
                                          "during")))

#just for checking!
#write.csv(dat_with_events, "Data/Processed_Data/dat_with_events.csv")




summary_by_event<- dat_with_events %>%
  group_by(event.id.y, phase, treat) %>%
  summarise(
    intensity = first(intensity),
    date = first(date),
    n.frames = n(),
    NumW_IQR = paste(unname(quantile(NumW, probs = c(0.25, 0.75))), collapse = "-"),
    NumW_min = min(NumW, na.rm = T), 
    NumW_max = max(NumW, na.rm = T),
    
    across(
      c(NumW, sd_head, Rank, File, 
        PVA, Mean_NN, mean_body_angle, Maj_axis),
      ~ paste0(
        round(mean(.x, na.rm = TRUE), 1),
        " (",
        round(sd(.x, na.rm = TRUE), 1),
        ")"
      )
    ),
    
    .groups = "drop"
  )


#sort neatly



summary_by_event_sorted <- summary_by_event %>%
  mutate(
    order_var = case_when(
      
      phase == "before" & treat == "more than 2 hours" ~ 1,
      phase == "before" & treat == "within 2hrs"       ~ 2,
      
      phase == "during"                                ~ 3,
      
      phase == "after" & treat == "within 2hrs"       ~ 4,
      phase == "after" & treat == "more than 2 hours" ~ 5,
      
      TRUE ~ 99
    )
  ) %>%
  arrange(event.id.y, order_var) %>%
  select(-order_var)





#~4d. Save-----
write.csv(summary_by_event, "Data/Processed_Data/obs_unfiltered_raw_summary.csv")



#5. Summary by type-----
summary_by_intensity<- dat_with_events %>%
  group_by(intensity, phase, treat) %>%
  summarise(
    intensity = first(intensity),
    date = first(date),
    n.frames = n(),
    NumW_IQR = paste(unname(quantile(NumW, probs = c(0.25, 0.75))), collapse = "-"),
    NumW_min = min(NumW, na.rm = T), 
    NumW_max = max(NumW, na.rm = T),
    
    across(
      c(NumW, sd_head, Rank, File, 
        PVA, Mean_NN, mean_body_angle, Maj_axis),
      ~ paste0(
        round(mean(.x, na.rm = TRUE), 1),
        " (",
        round(sd(.x, na.rm = TRUE), 1),
        ")"
      )
    ),
    
    .groups = "drop"
  )



