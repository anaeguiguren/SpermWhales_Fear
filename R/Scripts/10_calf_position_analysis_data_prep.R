#' ------------------------------------------------------------------------------
#' Title: Prepare data for calf position analysis by behavioral state
#' Author: Ana Eguiguren
#' Date: May 2026
#' Description: This script processes individual position metrics, assigns flight 
#'              session IDs, filters by behavioral events, and computes block-level
#'              summaries for GLM analysis
#' ------------------------------------------------------------------------------

# Clear workspace
rm(list = ls())

# 0. SETUP ----
# Load libraries
library(dplyr)
library(patchwork)
library(lubridate)
library(ggplot2)
library(tidyr)
library(here) 

#Load data and set parameters
load(here("R/labs_formations.RData"))
set.seed(12345)

LOCAL_TZ <- "Pacific/Galapagos"


# 1. FORMAT AND CLEAN BASE DATA ----
#~~1a. Convert timestamp to datetime ----
clust.whales <- clust.whales %>%
  mutate(frame_local_time_f = as.POSIXct(
    frame_local_time, 
    format = "%Y-%m-%d %H:%M:%S", 
    tz = LOCAL_TZ
  ))

#~~1b. Standardize age classes ----
# Convert NA to "Adult" and create coarse age classification
clust.whales <- clust.whales %>%
  mutate(
    age_class = if_else(is.na(age_class), "Adult", age_class),
    age_class_coarse = if_else(
      age_class %in% c("Adult", "MatMal", "Juv"), 
      "Adult", 
      "Calf"
    )
  )

# Verify age class distribution
cat("\nAge class distribution:\n")
print(table(clust.whales$age_class_coarse))

#~~1c. Sort by time and cluster ----
clust.whales <- clust.whales %>%
  arrange(frame_local_time_f, clust)


# 2. COMPUTE STANDARDIZED POSITIONS AND RANKS ----
#~~2a. Create unique cluster identifiers ----
clust.whales <- clust.whales %>%
  mutate(
    cluster_ID = paste(frame_ID, clust, sep = "_"), 
    f_clust_ID = as.numeric(factor(cluster_ID))
  )

#~~2b. Compute within-cluster standardized metrics ----
# Standardize positions relative to cluster min/max
# Keep only clusters with 4+ whales and at least one calf
clust.whales <- clust.whales %>%
  group_by(cluster_ID) %>%
  mutate(
    nWhal = n(),
    
    filep_f_std = 2 * ((filep_f - min(filep_f)) /
                         (max(filep_f) - min(filep_f))) - 1,
    
    # distance from center rank
    rankP_dist = abs(
      rank(rankP_l, ties.method = "average") - (nWhal + 1)/2
    ) - 0.5,
    
    # min-max scale to [0,1]
    rankP_min_std = (
      rankP_dist - min(rankP_dist)
    ) / (
      max(rankP_dist) - min(rankP_dist)
    ),
    
    # target sum if mean should equal 0.5
    target_sum = 0.5 * nWhal,
    
    # observed sum
    obs_sum = sum(rankP_min_std),
    
    # correction factor
    conv = target_sum / obs_sum,
    
    # corrected metric
    rankP_min_std = rankP_min_std * conv
  ) %>%
  filter(
    n() > 3,
    any(age_class_coarse == "Calf")
  ) %>%
  select(
    -rankP_dist,
    -target_sum,
    -obs_sum,
    -conv
  ) %>%
  ungroup()
cat("\nClusters after filtering: ", n_distinct(clust.whales$cluster_ID), "\n")


# 3. COMPUTE CLUSTER EDGE METRICS ----
#~~3a. Assign binary edge position ----
# Define "outer 25%" as those in top quartile of distance from center
clust.whales <- clust.whales %>%
  group_by(cluster_ID) %>%
  mutate(
    outer_25 = dist_cent >= quantile(dist_cent, probs = 0.75),
    dist_std = dist_cent / max(dist_cent),
    dist_scaled = (dist_cent - min(dist_cent)) / (max(dist_cent) - min(dist_cent))
  ) %>%
  ungroup()


# 4. MERGE BEHAVIORAL EVENT METADATA ----
#~~4a. Load and process event descriptions ----
event_desc <- read.csv("Input_Data/Event_Descriptions.csv", header = TRUE) %>%
  select(event.code, Intensity, species) %>%
  mutate(event.code = tolower(event.code))

#~~4b. Reclassify intensity levels by species ----
event_desc <- event_desc %>%
  mutate(
    intensity_sp = case_when(
      Intensity == "Neutral" ~ "Neutral",
      species == "Gm" ~ "Aversive_Gm",
      species == "Oo" ~ "Aversive_Oo",
      TRUE ~ NA_character_
    ),
    # Standardize intermediate to aversive
    Intensity = if_else(Intensity == "Intermediate", "Aversive", Intensity)
  )

#~~4c. Standardize event codes with underscores ----
event_desc <- event_desc %>%
  mutate(event.code = if_else(
    grepl("_", event.code),
    event.code,
    gsub(" ", "_", event.code)
  ))

#~~4d. Join event metadata to position data ----
clust.whales <- clust.whales %>%
  mutate(event.code = tolower(event)) %>%
  left_join(event_desc, by = "event.code")


# 5. FILTER OUT BUFFER PERIODS AROUND EVENTS ----
#~~5a. Load and process event timing data ----
event_times <- read.csv("Data/Processed_Data/Boris_Event_Times.csv", header = TRUE)

#~~5b. Get start/end times for each event ----
event_windows <- event_times %>%
  group_by(event) %>%
  summarize(
    start_date_time = first(start_date_time), 
    end_date_time = last(end_date_time), 
    duration_focal = sum(duration),
    .groups = "drop"
  ) %>%
  mutate(
    start_date_time = as.POSIXct(start_date_time, format = "%Y-%m-%d %H:%M:%S"),
    end_date_time = as.POSIXct(end_date_time, format = "%Y-%m-%d %H:%M:%S"),
    event.code = tolower(event)
  )

#~~5c. Join event metadata ----
event_windows <- event_windows %>%
  left_join(event_desc, by = "event.code")

#~~5d. Create 10-minute buffer zones ----
# Exclude observations within 10 min before start and 10 min after end
event_windows <- event_windows %>%
  mutate(
    start_buffer = start_date_time - minutes(10),
    end_buffer = end_date_time + minutes(10)
  )

#~~5e. Filter out observations in buffer zones ----
# More efficient approach: create exclusion intervals
cat("\nObservations before buffer filtering: ", nrow(clust.whales), "\n")

# Create a function to check if time falls in any exclusion zone
is_in_buffer <- function(time, windows) {
  any(
    (time >= windows$start_buffer & time <= windows$start_date_time) |
    (time >= windows$end_date_time & time <= windows$end_buffer)
  )
}

# Apply filtering (vectorized approach)
clust.whales <- clust.whales %>%
  filter(!sapply(frame_local_time_f, is_in_buffer, windows = event_windows))

cat("Observations after buffer filtering: ", nrow(clust.whales), "\n")


# 6. ASSIGN FLIGHT SESSION IDs ----
#~~6a. Load and process drone video coding data ----
drone_sessions <- read.csv("Data/Gal2023_Drone_Video_Coding_Social_Touch.csv") %>%
  filter(!is.na(Index)) %>%
  mutate(
    LocalTime_f = as.POSIXct(LocalTime, format = "%Y-%m-%d %H:%M", tz = LOCAL_TZ),
    flight_session_ID = paste0(Date, "_", session_id)
  ) %>%
  arrange(LocalTime_f)

#~~6b. Get start/end times for each flight session ----
session_summary <- drone_sessions %>%
  group_by(flight_session_ID) %>%
  summarize(
    session_start = first(LocalTime_f),
    session_end = last(LocalTime_f) + minutes(12),  # Max video length
    .groups = "drop"
  )

#~~6c. Assign session IDs to frames via time overlap ----
clust.whales <- clust.whales %>%
  left_join(
    session_summary, 
    by = join_by(
      frame_local_time_f >= session_start, 
      frame_local_time_f <= session_end
    )
  )

cat("\nObservations with assigned flight sessions: ", 
    sum(!is.na(clust.whales$flight_session_ID)), "\n")


# 7. CREATE BLOCKING VARIABLE BASED ON FLIGHT SESSION AND PHASE----
clust.whales<- clust.whales %>%
  arrange(flight_session_ID, frame_local_time_f)%>%
  group_by(flight_session_ID)%>%
  
  mutate(
    interaction = Intensity == "Aversive", # distinguish during vs not
    
    first_interaction = if (any(interaction, na.rm = TRUE)) {
      min(which(interaction))
    } else {
      NA_integer_
    },
    
    last_interaction = if (any(interaction, na.rm = TRUE)) {
      max(which(interaction))
    } else {
      NA_integer_
    },
    
    phase = case_when(
      is.na(first_interaction) ~ "no_interaction",
      
      interaction ~ "during",
      
      row_number() < first_interaction ~ "before",
      
      row_number() > last_interaction ~ "after"
    )
    
  ) %>%
  ungroup() 


#~~7a. Summarize clusters by frame and behavioral state ----
cluster_summary <- clust.whales %>%
  group_by(frame_ID, cluster_ID) %>%
  summarize(
    frame_local_time = first(frame_local_time_f), 
    behav_state = first(behav_state),
    .groups = "drop"
  ) %>%
  arrange(frame_local_time)






# 8. PREPARE FINAL DATASET ----
#~~8a. Select relevant columns ----
position_data <- clust.whales %>%
  select(
    # Identifiers
    date, 
    frame_ID, 
    cluster_ID, 
    flight_session_ID,
    
    # Time and behavior
    frame_local_time = frame_local_time_f,
    interaction,
    event,
    phase,
    
    # Individual characteristics
    age_class_coarse,
    
    # Position metrics
    NND,
    filep_f_std, 
    rankP_min_std, 
    dist_std, 
    outer_25,
    
    # Event metadata
    intensity_sp
  )



#make blocking average:

position_data <- position_data %>%
  mutate(
    block_id = paste0(flight_session_ID, 
                      "_", 
                      phase)
  )

#check levels of blocks
levels(factor(position_data$block_id))


#~~8b. Compute block-level summaries by age class ----
# This is the final output mentioned: average values per block-id-age class
block_summaries <- position_data %>%
  group_by(block_id, age_class_coarse) %>%
  summarize(
    # Sample sizes
    n_observations = n(),
    n_clusters = n_distinct(cluster_ID),
    
    # Temporal info
    first_time = first(frame_local_time),
    last_time = last(frame_local_time),
    
    # Behavioral context
    intensity_sp = first(intensity_sp),
    flight_session_ID = first(flight_session_ID),
    
    # Position metrics - means
    mean_NND = mean(NND, na.rm = TRUE),
    mean_filep_std = mean(filep_f_std, na.rm = TRUE),
    mean_rankP_std = mean(rankP_min_std, na.rm = TRUE),
    mean_dist_std = mean(dist_std, na.rm = TRUE),
    
    # Position metrics - standard deviations
    sd_NND = sd(NND, na.rm = TRUE),
    sd_filep_std = sd(filep_f_std, na.rm = TRUE),
    sd_rankP_std = sd(rankP_min_std, na.rm = TRUE),
    sd_dist_std = sd(dist_std, na.rm = TRUE),
    
    # Proportion in outer 25%
    prop_outer_25 = sum(outer_25)/ n(), 
    
    .groups = "drop"
  )

#get state:

block_summaries <- block_summaries %>%
  mutate(state = 
           ifelse(is.na(intensity_sp), "not", "aversive"))


cat("\nFinal block summaries:\n")
cat("  Total blocks: ", n_distinct(block_summaries$block_id), "\n")
cat("  Adult observations: ", sum(block_summaries$n_observations[block_summaries$age_class_coarse == "Adult"]), "\n")
cat("  Calf observations: ", sum(block_summaries$n_observations[block_summaries$age_class_coarse == "Calf"]), "\n")


# 9. SAVE OUTPUTS ----
#~~9a. Save individual-level data ----
saveRDS(
  position_data, 
  file = "Data/Processed_Data/calf_position_individual_data.RData"
)

#~~9b. Save block-level summaries ----
saveRDS(
  block_summaries,
  "Data/Processed_Data/calf_position_block_summaries.RData"
)

