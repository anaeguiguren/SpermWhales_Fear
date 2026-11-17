#' ---
#' Title: Subset PCA-R outcomes for cluster formation metrics GLMMs
#' Author: Ana Eguiguren
#' Date: May 2026
#' Description: This  PCA-R outcomes and filters only frames observed in
#' association with events, during, within 2 hrs, before/after 2hrs of 
#' an event.  
#' within 2 hrs of event
#' ---



#0. LOAD LIBRARIES-----
rm(list = ls())

#list of files
packages <- c("dplyr", "tidyr", 
              "lubridate")

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

#1. LOAD MEDIAN-CENTERED DATA-----

dat <- readRDS("Data/Processed_Data/PCA-R_median_centered_full_data_form.rds")
d <- read.csv("Input_Data/Event_Descriptions.csv")

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


#~2c. subset observations during and close to events----
dat_med_cent_2hr_window <- dat %>%
  inner_join(
    event_windows,
    by = join_by(
      frame_local_time_f >= start_include,
      frame_local_time_f <= end_include
    )
  ) %>%
  # Compute distance to event
  mutate(
    time_to_event = case_when(
      
      # observation occurs during event
      frame_local_time_f >= start_date_time &
        frame_local_time_f <= end_date_time ~ 0,
      
      # before event
      frame_local_time_f < start_date_time ~
        as.numeric(difftime(start_date_time,
                            frame_local_time_f,
                            units = "secs")),
      
      # after event
      frame_local_time_f > end_date_time ~
        as.numeric(difftime(frame_local_time_f,
                            end_date_time,
                            units = "secs"))
    )
  ) %>%
  
  # keep closest event per observation
  group_by(frame_local_time_f) %>%
  slice_min(time_to_event, n = 1, with_ties = FALSE) %>%
  ungroup()



#~2d. filter out frames too close to event (< 10 mins)----

obs_filtered <- dat_med_cent_2hr_window %>%
  filter(
    !(
      frame_local_time_f >= start_exclude &
        frame_local_time_f <= end_exclude &
        (
          frame_local_time_f < start_date_time |
            frame_local_time_f > end_date_time
        )
    )
  )



# lable as before/during/after
obs_filtered <- obs_filtered %>%
  mutate(
    phase = case_when(
      frame_local_time_f >= start_date_time & frame_local_time_f <= end_date_time ~ "during",
      frame_local_time_f <  start_date_time ~ "before",
      frame_local_time_f >  end_date_time   ~ "after"
    )
  )


#~2e. keep only events with at least one before/during phase----

# get summary of events by phase
events <- obs_filtered %>%
  group_by(event.code, phase) %>%
  summarise(n = n(), 
            date = first(date.x)) %>%
  arrange(as.Date(date))

# identify useful events
evs <- events %>%
  group_by(event.code)%>% 
  filter(any(phase == "during") & n_distinct(phase)>1) %>%
  select(event.code)%>%
  distinct(event.code)

#keep only corresponding observations
obs_filtered<- obs_filtered %>%
  filter(event.code %in% evs$event.code)


#3. SAVE------
#saveRDS(obs_filtered, "Data/Processed_Data/obs_filtered_med_cent.rds")

#4. CREATE A DATASET WITH ALL ROWS BUT WITH TIME WINDOW INFO----

#~4a. join dat with event windows to get event times----
#for each frame.
dat_with_events <- dat %>%
  left_join(event_windows, by = character()) %>%
  # Keep only event windows on the same day as the frame
  filter(as.Date(frame_local_time_f) == as.Date(start_date_time)) %>%
  # Calculate time differences
  mutate(
    time_to_start = as.numeric(difftime(frame_local_time_f, start_date_time, units = "secs")),
    time_to_end   = as.numeric(difftime(frame_local_time_f, end_date_time, units = "secs")),
    abs_time_to_event = pmin(abs(time_to_start), abs(time_to_end), na.rm = TRUE)
  ) %>%
  # For each frame, keep only the closest event (if multiple events per day)
  group_by(frame_local_time_f) %>%
  slice_min(abs_time_to_event, n = 1, with_ties = FALSE) %>%
  ungroup()

#~4b. Assign interaction_window label----
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

#~4c. assign a phase (before/after)----

# label as before/during/after
dat_with_events <- dat_with_events %>%
  mutate(
    phase = case_when(
      frame_local_time_f >= start_date_time & frame_local_time_f <= end_date_time ~ "during",
      frame_local_time_f <  start_date_time ~ "before",
      frame_local_time_f >  end_date_time   ~ "after"
    )
  )


#~4d. keep only events with at least two phases including during----
dat_with_events <- dat_with_events %>%
  filter(event.code %in% evs$event.code)

#remove frames within 10 minutes
dat_with_events <- dat_with_events %>%
  filter(interaction_window != "within 10 mins" |
           is.na(interaction_window))

#~4e. cleanup variable names ------
dat_with_events <- dat_with_events %>%
  mutate(intensity_sp = intensity_sp, 
         intensity = Intensity,
         treat = interaction_window)


#~4f. inherit intensity sp by event-----

#~4f. make event id the key column
dat_with_events <- dat_with_events %>%
  group_by(event.code) %>%
  mutate(behav_state_int = ifelse(behav_state_int == "not_interacting", 
                                  NA, behav_state_int))%>%
  fill(behav_state_int, .direction = "downup")%>%
  ungroup()




#~4d. Save-----
saveRDS(dat_with_events, "Data/Processed_Data/obs_unfiltered_med_cent.rds")

