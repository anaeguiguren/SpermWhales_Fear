library(tidyr)
library(tidyverse)
library(ggplot2)
library(data.table)
library(stringr)



# 1. auxiliary functions-----



count_behavior_events <- function(start_time, end_time, point_events) {
  filtered <- point_events %>%
    filter(Time > start_time, Time <= end_time)
  
  counts <- table(factor(filtered$Behavior, levels = c("contact", "jaw", "thrash", "tail_slap")))
  as.list(counts)
}


count_blows_events <- function(subject_id, start_time, end_time, blows) {
  filtered <- blows %>%
    filter(Subject == subject_id, Time > start_time, Time <= end_time)
  
  # Just counting "blow" events
  nrow(filtered)
}


# 2. read BORIS outputs ------

FolderPath = "Input_Data/Boris_Outputs/"

read.boris <- function(FolderPath){
  #Boris file names
  
  obs <-list.files(FolderPath, full.names = T, pattern = ".csv")
  
  #initialize dataframes------
  
  scan_samples_all <- tibble()
  focal_follows_all <- tibble()
  #Loop through observations----
  
  for(i in seq_along(obs)){
    d <- read.csv(obs[i], header = T)
    print(d$Observation.id[1])
    
    if(nrow(d)==0){
      scan_samples <- tibble()
      focal_binned <- tibble()
    }else{
      # a. Get relevant tables----
      
      visual_effort <- d %>% filter(Behavioral.category == "Sampling") # start and end of visual effort
      
      scan_samples <- d %>% filter(Behavior == "scan_behavs") # scan samples 
      
      point_events <- d %>% filter(Behavior == "jaw"| Behavior =="thrash" | Behavior =="contact" | Behavior == "tail_slap") #events
      
      focal_follows <- d %>% filter(Behavior == "focal_follow") # focal follows
      
      blows <- d %>% filter(Behavior == "blow")
      
      
      # b. clean visual effort-----
      #start and end of acceptable visual contact (excluding zoom in periods and periods with bad visibility)
      
      
      visual_effort <- visual_effort %>% 
        select(Event_ID = Observation.id , Behavior, n = Modifier..1, Behavior.type, Start.time = Start..s., Stop.time = Stop..s., Comments = Comment.start)%>%
        mutate(ff_id = paste(Event_ID, n, sep = "_"), n = as.numeric(n), Comments = as.character(Comments))
      
      visual_effort <-visual_effort %>%
        pivot_longer(
          cols = c(Start.time, Stop.time),
          names_to = "type", 
          names_pattern = "(Start|Stop)\\.time",
          values_to = "Time"
        )%>%
        arrange(Time)
      
      
      
      
      
      
      # c. clean scan samples-----
      
      
      
      scan_samples <-scan_samples %>% select(
        Event_ID = Observation.id, 
        Behavior, 
        Behavior.type,
        Time = Start..s., 
        N_sperms = Modifier..1, 
        N_other_sp = Modifier..2,
        N_nb = Modifier..3,
        N_c = Modifier..4,
        N_m = Modifier..5,
        Comments = Comment.start
      )%>%
        mutate(across(starts_with("N_"), as.character)) %>%
        mutate(across(starts_with("N_"), ~ if_else(.x == "None", "0", .x))) %>% #none are 0s
        mutate(across(starts_with("N_"), as.numeric)) %>%
        mutate(Comments = as.character(Comments))# 
      
      
      scan_samples<-bind_rows(visual_effort, scan_samples)%>%
        arrange(
          Time)
      
      
      
      # c. clean behaviours-----
      point_events <- point_events %>%
        select(
          Event_ID = Observation.id, 
          Behavior, Behavior.type, 
          Time = Start..s., 
          Comments = Comment.start
        )
      
      
      
      
      # d. count number of behaviours in time interval ----
      scan_samples<- scan_samples %>%
        mutate(Prev_Time = lag(Time, default = 0), 
               Time_interval = Time - Prev_Time)
      
      scan_samples$Time_interval[1]=NA
      
      
      #function for counting events:
      
      
      
      
      scan_samples <- scan_samples %>%
        rowwise() %>%
        mutate(counts = list(count_behavior_events(Prev_Time, Time, point_events))) %>%
        unnest_wider(counts, names_sep = "_") %>%
        ungroup()
      
      
      
      
      #e. estimate blow rates -----
      # clean:
      focal_follows <- focal_follows %>%
        select(
          Event_ID = Observation.id, 
          Behavior, Behavioral.category, 
          Subject,
          Start.Time = Start..s., 
          Stop.Time = Stop..s., 
          Comments = Comment.start
        )%>%
        mutate(Comments = as.character(Comments))
      
      blows <- blows %>%
        select(
          Event_ID = Observation.id, 
          Behavior, Behavioral.category, 
          Subject, 
          Time = Start..s., 
          Comments = Comment.start
        )
      
      
      
      
      #count for each individual how many blows were counted each minute:
      
      focal_binned <- focal_follows %>%
        rowwise() %>%
        mutate(
          #makes list of 60 sec bins:
          Time = list(seq(from = Start.Time, to = Stop.Time-1, by = 60 ))
        )%>%
        unnest(
          Time
        ) %>%
        arrange(Time)%>%
        select(-Start.Time, -Stop.Time)
      
      
      #get counts:
      focal_binned<- focal_binned %>%
        group_by(Subject)%>%
        mutate(Prev_Time = lag(Time, default = 0), 
               Time_interval = Time - Prev_Time)%>%ungroup()
      
      
      
      #function for counting events:
    
      
      #apply function
      
      focal_binned<-focal_binned %>%
        rowwise() %>%
        mutate(
          blows_in_period = count_blows_events(Subject, Prev_Time, Time, blows)
        ) %>%
        ungroup()%>%
        filter(Prev_Time!=0)
    }
    
    
    # f. join all the loops
    
    scan_samples_all<-bind_rows(scan_samples_all, scan_samples)
    focal_follows_all<-bind_rows(focal_follows_all, focal_binned)
    
  }
  
  #remove junk
  rm(blows, d, focal_binned, focal_follows, point_events, scan_samples, visual_effort)
  return(list(scan_samples = scan_samples_all, focal_follows = focal_follows_all))
  
}

# 3. estimate behavioural rates----
behav.rates <- function(scan_sample_data){
 
    #estimate general rates
    x<-scan_sample_data %>%
    mutate(jaw_rate = ifelse(Behavior.type =="STATE", NA, counts_jaw/N_sperms*(Time_interval/60)), 
           thrash_rate = ifelse(Behavior.type =="STATE",NA, counts_thrash/N_sperms*(Time_interval/60)), 
           tail_slap_rate = ifelse(Behavior.type =="STATE",NA,counts_tail_slap/N_sperms*(Time_interval/60)))
  
   #identify and flag zoomed in scan samples:
    zoom_periods <- x %>% 
      filter(Behavior == "zoom_in_state")%>%
      select(Event_ID, Time, type, Comments) %>%
      group_by(Event_ID)%>%
      mutate(pair_id = paste0(Event_ID, "_", cumsum(type== "Start")))%>%
      pivot_wider(names_from = type, values_from = Time) %>%
      filter(!is.na(Start) & !is.na(Stop)) %>%
      select(Event_ID, Start, Stop)
    
    
    #join to scan samples to flag bad ones:
   x %>%
      rowwise() %>%
      mutate(
        bad_samps = any(
          zoom_periods$Event_ID == Event_ID &
            Time >= zoom_periods$Start &
            Time < zoom_periods$Stop
        )
      )
    
  
}


#4. Read Sleap Labels -----
#read sleap labels raw
read.sleap <- function(txt.file){
  
  
  #txt.file <- "Sleap_Labels/Gal2023_Sleap_Project_raw_outputs_interspecies.txt"
  a <- read.delim(txt.file, header = F)
  
  #each row has all the instances for one frame and one video
  #for each frame, split text into one row/instance (whale)
  
  d <- tibble()
  
  for(i in seq_along(a$V1)){
    
    fr <-unlist(str_split(a[i,1], "Instance"))# all the instances for frame 1
    d.temp <- data.frame(
      filename = str_match(fr, "filename=\\s*(.*?)\\s*, shape")[,2],
      frame_idx = str_match(fr, "frame_idx=\\s*(.*?)\\s*, points")[,2],
      snout_pos = str_match(fr, "snout:\\s*(.*?)\\s*, dorsal")[,2],
      df_pos = str_match(fr, "dorsal:\\s*(.*?)\\s*, fluke")[,2],
      fluke_pos = str_match(fr, "fluke:\\s*(.*?)\\s*], track")[,2]
    )
    
    d.temp <- d.temp[complete.cases(d.temp),]
    d<-rbind(d, d.temp)
  }
  
  
  
  
  # Sort the original dataframe based on video and frame index 
  d <-d %>%
    arrange(filename, as.numeric(frame_idx))

  # Print the sorted dataframe
  
  
  
  rm(a, d.temp, fr)
  return (d)
  
  
  
}


# 5. Join labels to drone altitude ------
# match to drone data for analyses
#logs.file = "C:/Users/balae/Documents/Galapagos_Drone_Logs/Gal2023_Drone_Flight_Logs.csv"
#labels = labs



#data = labs
# 2. getSRTAltitude -----
getSrtAltitude <- function(data){
  #read in all flight srt data
  require(stringr)
  drone_srt_files <- read.csv("Data/Drone_Logs/Gal2023_Drone_Flight_Logs_srt.csv", header = T)
  
  #make matching file name
  drone_srt_files <- drone_srt_files %>%
    mutate(mp4_file = basename(videoFile))
  
  #obtain snapshot file type: 
  data <- data %>% 
    mutate(  # identify vlc and other file snapshot types
      mp4_file = basename(video),
      ss_sec = floor(frame_idx/29.97)#snapshot time in seconds, rounded down
    )
  
  #join based on file name and time in video
  
  
  
  # Join data with drone_srt_log to get altitude based on mp4_file & time match
  data <- data %>%
    left_join(drone_srt_files, 
              by = c("mp4_file" = "mp4_file", "ss_sec" = "videoTime")) %>%
    select(date, group_ID, Q, behav_state, frame_local_time, interaction, event,
           video, frame_idx, instance_id, 
           x_snout, y_snout,
           x_dorsal, y_dorsal, 
           x_fluke, y_fluke,
           snapshot_path, ss_sec, imageWidth, droneDateTime, droneDate, 
           droneFlightLatitude, droneFlightLongitude, droneAltitude, 
           age_class, whale_length
           )
  

  # View result
  return(data)
  
}


#5. Correct altitude with correction factor and launch height----
altitudeASL <- function(boat_height =  1.03 - 0.24, 
                        launch.chest = 1.4, 
                        camera.height = 0.045, 
                        altitude.raw){
  altitude.fix = altitude.raw+launch.chest+camera.height+boat_height
  altitude.c = 1.40 + altitude.fix*1.017
  return(altitude.c)
  
}



#8. add.event.info------

add.event.info <- function(events, data){
  events_expanded <- events %>%
    mutate(mp4_files =  str_extract_all(Media.files.Live.observation, "[^/\\\\]+\\.MP4"))%>%
    unnest(mp4_files)
  
  data %>% 
    left_join(events_expanded, by = c("mp4_file" = "mp4_files"))
  
}
