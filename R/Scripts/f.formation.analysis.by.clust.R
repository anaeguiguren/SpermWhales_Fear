#Formation analysis:
f.formation.analysis.clust <- function(data){
  
  #' @title formation.analysis
  #' @description reads labeled frames with whale positions and outputs formation analysis statistics
  #' @param labeled_frame_path full path to where the labeled frames data are stored (character)
  
  
  require(circular);require(dplyr);require(amap); require(dbscan)
  # source("formation_analysis/f.clumping.R")
  
  # Data ----
  frames <- levels(as.factor(data$frame_ID))
  
  #initialize empty data frames-----
  # for frames, by cluster
  aa <- data.frame(date = as.character(),
                   group_ID = as.character(),
                   Q = as.character(), 
                   behav_state = as.character(),
                   frame_local_time = as.character(), 
                   interaction = as.character(), 
                   event = as.character(),
                   filename = as.character(),
                   time_stamp_sec = as.character(),
                   frame_ID = as.character(),
                   local_time_stamp = as.numeric(),
                   droneH = as.numeric(),
                   conv = as.numeric(),
                   
                   #cluster level things
                   clust = as.numeric(), #cluster number
                   NumW = as.numeric(), #number of whales in cluster
                   mean_head = as.numeric(),# mean heading in degrees
                   sd_head = as.numeric(), # standard deviation of heading in degrees
                   rho_head = as.numeric(), #rho circ stat
                   mean_body_angle = as.numeric(), #mean body angle at dorsal fin 
                   Rank = as.numeric(), #rank length in meters (width)
                   File = as.numeric(), # File length in meters 
                   File_Rank_ratio = as.numeric(), #File rank ratio
                   Maj_axis = as.numeric(), # angle of PC1 with respect to average heading
                   PVA = as.numeric(), # proportion of variance explained by principal component 1
                   Mean_NN = as.numeric()
  )
  #for whales:
  bb <- data.frame(
    date = as.character(),
    group_ID = as.character(),
    Q = as.character(), 
    behav_state = as.character(),
    frame_local_time = as.character(), 
    interaction = as.character(), 
    event = as.character(), 
    filename =as.character(), 
    time_stamp_sec = as.character(),
    
    FrameID = as.character(), #Frame ID
    local_time_stamp = as.numeric(),
    instance_id = as.character(), #stable instance ID
    clust = as.numeric(), #cluster id within frame
    nw = as.character(),# whale number
    length = as.numeric(), # length in meters
    NND = as.numeric(), #Distance to nearest neighbor
    rankF = as.numeric(), # position along rank in meters
    fileD = as.numeric(), # position along file in meters
    headSDf = as.numeric(), # individual heading
    headV = as.numeric(), # deviation from mean heading in degrees
    age_class = as.character(),#manually assigned age class
    whale_length = as.numeric(),# 'good' length measurement
    body_angle = as.numeric(),
    dist_cent = as.numeric()
  )
  # loop through each frame ----
  for(f in seq_along(frames)){
    
    # stuff needed for further calculations
    frame_f <- data %>%
      filter(frame_ID == frames[f])
    
    # drone altitude:
    droneH = frame_f[1, "c_altitude"]
    
    image.width = frame_f[1, "imageWidth"]
    
    # Conversion factor (GSD * altitude)
    
    alpha = ifelse(image.width == 3840, yes = 0.000328, no = 
                     ifelse(image.width == 1920, yes = 0.000656, no = NA))
    
    conv = alpha * droneH 
    
    
    # substracts dorsal fin coordinates - head coordinates to generate a triangle from which to calculate hypothenuse direction
    azd = data.frame(
      cbind(
        frame_f$x_dorsal-frame_f$x_snout,
        frame_f$y_dorsal-frame_f$y_snout
      )
    )# subtracts  fluke coordinates - head coordinates
    
    headings = atan2(-azd[,1],azd[,2]) #estimate headings (in radians)
  

    # substracts fluke  coordinates - dorsal fin coordinates to generate a triangle from which to calculate hypothenuse  
    #and add for length
    
    
    azd.f = data.frame(
      cbind(
        frame_f$x_fluke-frame_f$x_dorsal,
        frame_f$y_fluke-frame_f$y_dorsal
      )
    )# subtracts  fluke coordinates - head coordinates
    
    lengths = (sqrt(azd.f[,1]^2 + azd.f[,2]^2)+ sqrt(azd[,1]^2 + azd[,2]^2)) * conv # length in meters from fluke to head
    
    # whale body angles 
    
    whale_vectors <- frame_f %>%
      mutate(ux = x_snout - x_dorsal, #snout to dorsal line
             uy = y_snout - y_dorsal, 
             vx = x_fluke - x_dorsal, #fluke to dorsal line
             vy = y_fluke - y_dorsal, 
             angle_rad = acos(
               (ux * vx + uy * vy)/ #dot product
                 (sqrt(ux^2 + uy^2) * sqrt(vx^2 +vy^2)) #magnitudes
             ),
             angle_deg = 180 - angle_rad *180/pi
             )
    
    body_angles <- as.circular(whale_vectors$angle_deg, units = "degrees")
    
    
    # Orientation Analysis:
    
    # Rotate whales:
    whpos = frame_f %>% select(x_snout, y_snout) 
    
    whpos.mean = apply(whpos, 2, function(x) mean(x,na.rm = T)) #mean head position
    
    whpos = cbind(whpos[,1] -whpos.mean[1] , whpos[,2] -whpos.mean[2] ) * conv #center positions around mean and multiply by conv
    
    
    
    # by cluster
    
    c <- dbscan(whpos, eps = 30, minPts = 1) # assign cluster based on 30 m distance
    whpos<- data.frame(cbind(whpos, cluster = c$cluster)) #add cluster info
    
    
    # center by cluster:
    
    whpos <- whpos %>%
      group_by(cluster) %>%
      mutate(mean_x = mean(V1), 
             mean_y = mean(V2))%>%
      ungroup()%>%
      mutate(
        x = V1 - mean_x, 
        y = V2 - mean_y
      )
    
    
    #convert to polar coordinates:
    rho = sqrt(whpos$x^2 + whpos$y^2)# vector length in meters
    th = atan2(- whpos$y, whpos$x) # angle in radians
    
    
    whpos <- as.data.frame(cbind(whpos, headings))
    
    whs <- data.frame(x = whpos$x, 
                      y = whpos$y, 
                      cluster = whpos$cluster, 
                      headings = headings, 
                      lengths = lengths,
                      body_angles = body_angles,
                      rho, 
                      th)#initialize nearest neighbor
    
    whs<- whs %>%
      group_by(cluster) %>%
      mutate(
        mean_heading_c = mean(as.circular(headings, units ="radians"), na.rm = T),#mean heading by cluster
        xu = rho * as.numeric(cos(th + mean_heading_c)), # rotated lengths
        yu = rho * as.numeric(sin(th + mean_heading_c)), 
        rankP_l = rank(xu, ties.method = "min"), #rank left
        rankP_r = rank(-xu, ties.method = "min"), #rank position right
        rankP_min = pmin(rankP_l, rankP_r), #minimum rank position
        
        fileP_b = rank(yu, ties.method = "min"), #file position (back)
        fileP_f = rank(-yu, ties.method = "min"), #file position (lead)
        dist_cent = sqrt(xu^2 + yu^2)
        
      )%>%
      ungroup()%>%
      mutate(NN = NA)#initialize position within cluster
    
    #initialize by cluster
    # Number of whales
    
    cluster_info<- whs %>%
      group_by(cluster) %>%
      summarize(NumW = n(), 
                mean_heading= as.numeric(deg(mean(as.circular(headings, types = "angles", units ="radians"), na.rm = T))),
                sd_heading = as.numeric(deg(sd.circular(headings, types = "angles", units ="radians", na.rm = T))),
                rho_circular = as.numeric(rho.circular(as.circular(headings, units = "radians"))),
                body_angle = as.numeric(mean(body_angles, na.rm = T)),
                Rank = max(xu, na.rm = T) - min(xu, na.rm = T), #calculate width of rank (in meters) perpendicular to mean heading
                File = max(yu, na.rm = T) - min(yu, na.rm = T), #calculate length of file (in meters) along to mean heading
                File_Rank_ratio = NA,
                PVA = NA, #initialize pva , mean NN, and Maj_axis, Angle
                Mean_NN = NA, 
                Maj_axis = NA) 
    
    
    for(i in seq_along(cluster_info$cluster)) {
      pca.dat <- whs %>% filter(cluster == i) %>% select(x = xu, y = yu)
      pca.out <- pca(pca.dat, reduce = F, center = F)
      
      
      
      pc1 <- pca.out$loadings[, "Comp 1"]
      
      # enforce orientation (in case major axis is flipped)
      if (pc1["x"] < 0) pc1 <- -pc1
      
      
      if(length(coef) >2){
        coeff <- -coeff
      }
      
      
      cluster_info$Maj_axis[which(cluster_info$cluster==i)] = acos(
        abs(sum(pc1 * c(x = 0, y = 1)))
      ) * 180 / pi
      
      
      
      score = pca.out$scores # scores = positions of each point
      
      latent = pca.out$sdev^2 # latent = principal component variances/eigenvalues
      
      cluster_info$PVA[which(cluster_info$cluster==i)] = latent[1]/sum(latent) #Proportion of Variance accounted for by major axis (PC1)
      
      
      cluster_info$File_Rank_ratio[which(cluster_info$cluster==i)] = cluster_info$File[which(cluster_info$cluster==i)]/
        cluster_info$Rank[which(cluster_info$cluster==i)]
      
      #nearest neigbor distances:
      
      xxu = matrix(rep(pca.dat$x, length(pca.dat$x)), nrow = length(pca.dat$x), ncol=length(pca.dat$x), byrow=F)
      yyu = matrix(rep(pca.dat$y, length(pca.dat$y)), nrow = length(pca.dat$y), ncol=length(pca.dat$y), byrow=F)
      
      #Distances between whales in meters:
      ddu = sqrt(((xxu-t(xxu))^2) + ((yyu-t(yyu))^2)) 
      diag(ddu) = NA
      
      #Mean nearest neighbour distance:
      #average per cluster
      cluster_info$Mean_NN[which(cluster_info$cluster==i)] = mean(apply(ddu, 1, function(x) min(x, na.rm = T)), na.rm = T)
      
      #per whale within cluster
      whs$NN[which(whs$cluster == i)] = apply(ddu, 1, function(x) min(x, na.rm = T))
      
      
    }
    
    
    
    #Clump = clarke-evans clumpling factor in a 20 x 20 m window
    
    #whl_pp <- ppp(whs$xu, whs$yu, window = owin(c(-10, 10), c(-10, 10)))
    #ce <-clarkevans(whl_pp, correction = "Donnelly") 
    
    
    
    #initialize tmp. dataframe by frame
    aa.temp <- tibble(
      date = as.character(),
      group_ID = as.character(),
      Q = as.character(), 
      behav_state = as.character(),
      frame_local_time = as.character(), 
      interaction = as.character(), 
      event = as.character(),
      filename = as.character(),
      time_stamp_sec = as.character(),
      frame_ID = as.character(),
      local_time_stamp = as.character(),
      droneH = as.numeric(),
      conv = as.numeric(),
      
      #cluster level tings
      clust = as.numeric(), #cluster number
      NumW = as.numeric(), #whales in cluster
      mean_head  = as.numeric(),# mean heading in degrees
      sd_head = as.numeric(), # standard deviation of heading in degrees
      rho_head = as.numeric(),
      mean_body_angle = as.numeric(),
      Rank = as.numeric(), #rank length in meters (width)
      File = as.numeric(), # File length in meters 
      File_Rank_ratio = as.numeric(), # File to Rank ratio (long vs wide)
      Maj_axis = as.numeric(), # angle of PC1 with respect to average heading
      PVA = as.numeric(), 
      Mean_NN = as.numeric()
    )
    
    
    aa.temp.temp <- tibble(
      date = frame_f[1, "date"],
      group_ID = frame_f[1, "group_ID"],
      Q = frame_f[1, "Q"], 
      behav_state = frame_f[1, "behav_state"],
      frame_local_time = frame_f[1, "frame_local_time"], 
      interaction = frame_f[1, "interaction"], 
      event = frame_f[1, "event"],
      filename = frame_f[1, "video"], 
      time_stamp_sec = frame_f[1, "ss_sec"],
      frame_ID = frame_f[1, "frame_ID"],
      droneH,
      conv,
      
      #cluster level tings
      clust = cluster_info$cluster, #cluster id
      NumW = cluster_info$NumW, #whales in cluster
      mean_head = cluster_info$mean_heading,# mean heading in degrees
      sd_head = cluster_info$sd_heading, # standard deviation of heading in degrees, 
      rho_head = cluster_info$rho_circular,
      mean_body_angle = cluster_info$body_angle,
      Rank = cluster_info$Rank, #rank length in meters (width)
      File = cluster_info$File, # File length in meters 
      File_Rank_ratio = cluster_info$File_Rank_ratio, # file rank ration
      Maj_axis = cluster_info$Maj_axis, # angle of PC1 with respect to average heading
      PVA = cluster_info$PVA, 
      Mean_NN = cluster_info$Mean_NN
      
    )
    
    aa.temp <- rbind(aa.temp, aa.temp.temp)
    
    #for each whale
    
    bb.temp <- tibble(
      date = as.character(),
      group_ID = as.character(),
      Q = as.character(), 
      behav_state = as.character(),
      frame_local_time = as.character(), 
      interaction = as.character(), 
      instance_id = as.character(), #stable instance ID
      event = as.character(),
      filename = as.character(),
      time_stamp_sec = as.character(),
      frame_ID = as.character(),
      local_time_stamp = as.character(),
      
      #cluster stuff
      clust = as.numeric(), #cluster ID
      nw=as.numeric(), #whale number
      length = as.numeric(), # whale length (in meters)
      NND = as.numeric(), #Nearest neighbour distance
      rankD = as.numeric(),# Rank displacement (in meters)
      fileD= as.numeric(), # File displacement (in meters)
      headV = as.numeric(),#heading variation from the mean in degrees 
      age_class = as.character(),#manually assigned age class
      whale_length = as.numeric(),# 'good' length measurement
      rankP_l = as.numeric(), 
      rankP_r = as.numeric(), 
      rankP_min = as.numeric(), 
      fileP_b = as.numeric(), 
      fileP_f = as.numeric(), 
      body_angle = as.numeric(),
      dist_cent = as.numeric()
      
      
      
      
    )
    #for each whale in the frame:
    
    bb.temp.temp <- tibble(
      date = frame_f[1, "date"],
      group_ID = frame_f[1, "group_ID"],
      Q = frame_f[1, "Q"], 
      behav_state = frame_f[1, "behav_state"],
      frame_local_time = frame_f[1, "frame_local_time"], 
      interaction = frame_f[1, "interaction"], 
      event = frame_f[1, "event"],
      filename = frame_f[1, "video"], 
      time_stamp_sec = frame_f[1, "ss_sec"],
      frame_ID = frame_f[1, "frame_ID"],
      instance_id = frame_f[i, "instance_id"], #stable instance ID
      
      nw=i, #whale number
      clust = whs$cluster, #cluster ID
      length = whs$lengths, # whale length (in meters - measured locally)
      NND = whs$NN, #Nearest neighbour distance
      rankD = whs$xu,# Rank displacement (in meters)
      fileD= whs$yu, # File displacement (in meters)
      
      rankP_l = whs$rankP_l,# left flank pos
      rankP_r = whs$rankP_r,# right flank pos
      rankP_min = whs$rankP_min, #min rank position
        
      filep_f = whs$fileP_f,
      fileP_b = whs$fileP_b,
      
      
      headV = deg(whs$headings) - as.numeric(deg(whs$mean_heading_c)),#heading variation from the mean in degrees  
      age_class = frame_f[, "age_class"],#manually assigned age class
      whale_length = frame_f[, "whale_length"],# 'good' length measurement
      body_angle = whs$body_angles, 
      dist_cent = whs$dist_cent
    )
    
    bb.temp <- rbind(bb.temp, bb.temp.temp)
    
    
    
    
    
    
    aa <- rbind(aa, aa.temp)
    bb <- rbind(bb, bb.temp)
    
  }
  
  aa <- aa[order(aa$filename, as.numeric(aa$time_stamp_sec)), ]
  bb <- bb[order(bb$filename, as.numeric(bb$time_stamp_sec)), ]
  out <- list(aa, bb)
  names(out) <-c ("frames", "whales")
  
  return(out)
  
  
}

