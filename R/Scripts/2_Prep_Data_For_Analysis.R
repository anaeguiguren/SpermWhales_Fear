#' ---
#' Title: Prepare data for multivariate analysis
#' Author: Ana Eguiguren
#' Date: May 2026
#' Description: This script processes formation metrics and obtains
#' boat-based estimates of sociality index for downstream analysis
#' ---


#0. LOAD DATA AND LIBRARIES -----
rm(list = ls())

library(here)
library(dplyr)
library(purrr)
library(tidyr)
library(lubridate)

#source("R/Scripts/2_Prepare_Formation_Data.R")
load(here("R/labs_formations.RData"))

# 1. CLEAN AND FILTER FRAME-LEVEL DATA -----


# Build a compact dataset with the variables needed for analysis.
clust.frames.comp <- clust.frames %>%
  mutate(
    interaction_state = ifelse(is.na(interaction), "not", interaction),
    event_state       = ifelse(is.na(interaction), "not", event),
    event             = ifelse(is.na(event), "not", event)
  ) %>%
  dplyr::select(
    date,
    frame_local_time,
    frame_ID,
    group_ID,
    behav_state,
    event,
    interaction_state,
    event_state,
    NumW,
    mean_head,
    sd_head,
    Rank,
    File,
    PVA,
    Mean_NN,
    mean_body_angle,
    Maj_axis
  ) %>%
  filter(complete.cases(.))

# Keep only frames that are either outside an event or have 
# a documented interaction.
clust.frames.comp <- clust.frames.comp %>%
  filter(event == "not" | interaction_state != "not")


#2. ASSIGN INTERACTION INTENSITY LABELS----
#~2a. Load event metadata and classify interaction intensity----
d <- read.csv(here("Input_Data/Event_Descriptions.csv"), header = TRUE)

d <- d %>%
  dplyr::select(event.code, event.id, Intensity, species) %>%
  mutate(event.code = tolower(event.code)) %>%
  mutate(
    intensity_sp = ifelse(
      Intensity == "Neutral",
      "Neutral",
      ifelse(species == "Gm", "Aversive_Gm", "Aversive_Oo")
    )
  )

#~2b. Join intensity labels to the frame-level dataset----
clust.frames.comp <- clust.frames.comp %>%
  mutate(event.code = tolower(event))

dat <- left_join(clust.frames.comp, d, by = "event.code")

# Use the event-derived label when available.
dat <- dat %>%
  mutate(behav_state = ifelse(is.na(intensity_sp), behav_state, intensity_sp))

# Remove temporary objects no longer needed.
rm(clust.frames, clust.frames.comp, clust.whales, clust.form, labs)



#3. COMPUTE SOCIALITY INDEX FROM BOAT BASED DATA-----


#~3a. Load 5-minute cluster summaries and convert timestamps.
clu <- read.csv(here("Data/Processed_Data/summary_5min_behavioural_state.csv"))

LOCAL_TZ <- "Pacific/Galapagos"

clu <- clu %>%
  mutate(
    interval_start = as.POSIXct(interval_start, format = "%Y-%m-%d %H:%M:%S", tz = LOCAL_TZ),
    interval_end  = as.POSIXct(interval_end,   format = "%Y-%m-%d %H:%M:%S", tz = LOCAL_TZ)
  )
  

#~3b. Convert frame timestamps for window matching----
dat <- dat %>%
  mutate(
    frame_local_time_f = as.POSIXct(frame_local_time, 
                                  format = "%Y-%m-%d %H:%M:%S", tz = LOCAL_TZ)
  )





#~3c. Extract cluster statistics in a ±15-minute window around each frame----

WINDOW_SEC <- 15 * 60

dat <- dat %>%
  mutate(
    clu_stats = map(frame_local_time_f, function(t) {
      clu_window <- clu %>%
        filter(
          interval_start >= t - WINDOW_SEC,
          interval_start <= t + WINDOW_SEC
        )

      if (nrow(clu_window) == 0) {
        return(tibble(
          n_sampling_intervals = 0L,
          total_adults         = NA_real_,
          mean_cluster_size    = NA_real_,
          max_cluster_size     = NA_real_,
          total_fluke_ups      = NA_real_
        ))
      }

      tibble(
        n_sampling_intervals = n_distinct(clu_window$interval_start),
        total_adults         = sum(clu_window$Ad, na.rm = TRUE),
        mean_cluster_size    = mean(clu_window$Ad, na.rm = TRUE),
        max_cluster_size     = max(clu_window$Ad, na.rm = TRUE),
        total_fluke_ups      = sum(clu_window$fluke_ups, na.rm = TRUE)
      )
    })
  ) %>%
  unnest(clu_stats)


#~3d. compute fluke-up rate----
b_rate <- function(tally, whales, scan_dur) {
  tally / (whales * scan_dur)
}

dat <- dat %>%
  mutate(fluke.r = b_rate(tally = total_fluke_ups, whales = total_adults, scan_dur = n_sampling_intervals))

#~3e. Estimate the sociality index----
dat <- dat %>%
  mutate(soc = -0.240 + (0.06 * max_cluster_size) - (3.036 * fluke.r))


#remove non-important data sets

rm( clu)

#4. IDENTIFY FRAMES WITHIN FLIGHT SESSIONS -----

#~4a. get session id from video coding -----
dron <- read.csv(here("Data/Gal2023_Drone_Video_Coding_Social_Touch.csv"))

# format time stamp & name flight sessions uniquely:
dron <- dron %>%
  filter(!is.na(Index))%>%
  mutate(LocalTime_f = as.POSIXct(LocalTime, 
                                  format = "%Y-%m-%d %H:%M", tz = LOCAL_TZ), 
         flight_session_ID = paste0(Date, "_", session_id)) %>%
  arrange(LocalTime_f)


#get start end times for each flight session
session_summary <- dron %>%
  group_by(flight_session_ID) %>%
  summarize(session_start = first(LocalTime_f),
            session_end = last(LocalTime_f)+minutes(12))#maximum video length


#~4a. assing session id to frames----

dat <- dat %>%
  left_join(session_summary, 
            by = join_by(
              frame_local_time_f >= session_start, 
              frame_local_time_f <= session_end
            )
            )


#~4b. estimate rank/file ratio-----
dat <- dat %>%
  mutate(RFRatio = Rank/File)
rm(dron, session_summary)
