#Formation analysis:
data = labs
f.formation.analysis <- function(data){
  
  #' @title formation.analysis
  #' @description reads labeled frames with whale positions and outputs formation analysis statistics
  #' @param labeled_frame_path full path to where the labeled frames data are stored (character)

  
require(circular);require(dplyr);require(amap)
# source("formation_analysis/f.clumping.R")

  # Data ----
   frames <- levels(as.factor(data$frame_ID))
   
    #initialize empty data frames-----
      # for frames
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
                     NumW = as.numeric(),
                     mean_head = as.numeric(),# mean heading in degrees
                     sd_head = as.numeric(), # standard deviation of heading in degrees
                     Rank = as.numeric(), #rank length in meters (width)
                     File = as.numeric(), # File length in meters 
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
                    nw = as.character(),# whale number
                    length = as.numeric(), # length in meters
                    NND = as.numeric(), #Distance to nearest neighbor
                    rankF = as.numeric(), # position along rank in meters
                    fileD = as.numeric(), # position along file in meters
                    headSDf = as.numeric(), # individual heading
                    headV = as.numeric(), # deviation from mean heading in degrees
                    age_class = as.character(),#manually assigned age class
                    whale_length = as.numeric()# 'good' length measurement
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
          
          # Number of whales
          NumW = frame_f %>% 
            nrow()
            
      
          # substracts dorsal fin coordinates - head coordinates to generate a triangle from which to calculate hypothenuse direction
          azd = data.frame(
            cbind(
              frame_f$x_dorsal-frame_f$x_snout,
              frame_f$y_dorsal-frame_f$y_snout
              )
            )# subtracts  fluke coordinates - head coordinates
      
          headings = atan2(-azd[,1],azd[,2]) #estimate headings (in radians)
          

          # substracts fluke  coordinates - dorsal fin coordinates to generate a triangle from which to calculate hypothenuse  and add for length
    
      
          azd.f = data.frame(
            cbind(
              frame_f$x_fluke-frame_f$x_dorsal,
              frame_f$y_fluke-frame_f$y_dorsal
            )
          )# subtracts  fluke coordinates - head coordinates
          
          lengths = (sqrt(azd.f[,1]^2 + azd.f[,2]^2)+ sqrt(azd[,1]^2 + azd[,2]^2)) * conv # length in meters from fluke to head
          
          # Orientation Analysis:

                  # Rotate whales:
                  whpos = frame_f %>% select(x_snout, y_snout) 
                  
                  whpos.mean = apply(whpos, 2, function(x) mean(x,na.rm = T)) #mean head position
                  
                  whpos = cbind(whpos[,1] -whpos.mean[1] , whpos[,2] -whpos.mean[2] ) * conv #center positions around mean and multiply by conv
                  #convert to polar coordinates:
                  rho = sqrt(whpos[,1]^2 + whpos[,2]^2)# vector length in meters
                  th = atan2(-whpos[,2],whpos[,1]) # angle in radians
                  
                  #rotate whales so that mean heading is up ^, and then re-convert to carthesian coordinates:
                  xu = rho * cos(th + mean.circular(headings, types = "angles", units ="radians", na.rm = T))
                  yu = rho * sin(th + mean.circular(headings, types = "angles", units ="radians", na.rm = T))
                  
                  #calculate width of rank (in meters) perpendicular to mean heading:
                  Rank = (max(xu, na.rm = T)-min(xu, na.rm = T))
                  
                  #calculate length of file (in meters) along to mean heading:
                  File = (max(yu, na.rm = T)-min(yu, na.rm = T))
  
                  pca.dat <- data.frame(xu, yu); 
                  pca.dat<- pca.dat[complete.cases(pca.dat),]
                  
                  #PCA analysis & clumping factor
                  #estimate pca coefficients manually because princomp 
                  #keeps doing silly things with signs
                  
                  ####### PCA #######
                  #pca.a <- princomp(pca.dat)
                  #pca <- prcomp(pca.dat)
                  
                  #coeff <- pca$rotation
                  
                  
                  
                  pca <- pca(pca.dat, reduce = F, center = F)
                  coeff <- pca$loadings
                  
                  
                  if(length(which(coeff < 0)) >2){
                    coeff <- -coeff
                  }
                  
                  ###############

                  Maj_axis = deg(atan2(coeff[1,2],coeff[1,1]))#major axis = orientation between the major axis and the whale heading (in degrees)
                  
                  
                  score = pca$scores # scores = positions of each point
                  
                  latent = pca$sdev^2 # latent = principal component variances/eigenvalues
                  

                  PVA = latent[1]/sum(latent) #Proportion of Variance accounted for by major axis (PC1)
                  
                  #calculate distances between whales:
                  #make x and y matrices:
                  xxu = matrix(rep(xu, length(xu)), nrow = length(xu), ncol=length(xu), byrow=F)
                  yyu = matrix(rep(yu, length(yu)), nrow = length(yu), ncol=length(yu), byrow=F)
                  
                  #Distances between whales in meters:
                  ddu = sqrt(((xxu-t(xxu))^2) + ((yyu-t(yyu))^2)) 
                  diag(ddu) = NA
                  
                  #Mean nearest neighbour distance:
                  Mean_NN = mean(apply(ddu, 1, function(x) min(x, na.rm = T)), na.rm = T)
                
                  
                  #Clump = clarke-evans clumpling factor in a 20 x 20 m window
                  
                  #whl_pp <- ppp(xu, yu, window = owin(c(-10, 10), c(-10, 10)))
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
                    NumW = as.numeric(),
                    mean_head  = as.numeric(),# mean heading in degrees
                    sd_head = as.numeric(), # standard deviation of heading in degrees
                    Rank = as.numeric(), #rank length in meters (width)
                    File = as.numeric(), # File length in meters 
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
                    NumW,
                    mean_head = as.numeric(deg(mean.circular(headings, types = "angles", units ="radians", na.rm = T))),# mean heading in degrees
                    sd_head = as.numeric(deg(sd.circular(headings, types = "angles", units ="radians", na.rm = T))), # standard deviation of heading in degrees
                    Rank = Rank, #rank length in meters (width)
                    File = File, # File length in meters 
                    Maj_axis = Maj_axis, # angle of PC1 with respect to average heading
                    PVA = unname(PVA), 
                    Mean_NN = Mean_NN
                    
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
                    nw=as.numeric(), #whale number
                    length = as.numeric(), # whale length (in meters)
                    NND = as.numeric(), #Nearest neighbour distance
                    rankD = as.numeric(),# Rank displacement (in meters)
                    fileD= as.numeric(), # File displacement (in meters)
                    headV = as.numeric(),#heading variation from the mean in degrees 
                    age_class = as.character(),#manually assigned age class
                    whale_length = as.numeric()# 'good' length measurement
                  )
                  #for each whale:
                  for(i in seq_along(rho)){

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
                      length = lengths[i], # whale length (in meters)
                      NND = min(ddu[i,], na.rm = T), #Nearest neighbour distance
                      rankD = xu[i]*conv,# Rank displacement (in meters)
                      fileD= yu[i]*conv, # File displacement (in meters)
                      headV = deg(headings[i]) - as.numeric(deg(mean.circular(headings, types = "angles", units ="radians", na.rm = T))),#heading variation from the mean in degrees  
                      age_class = frame_f[i, "age_class"],#manually assigned age class
                      whale_length = frame_f[i, "whale_length"]# 'good' length measurement
                      )
                    
                    bb.temp <- rbind(bb.temp, bb.temp.temp)
                    
                    
                  
                  
                  
                }
          
      aa <- rbind(aa, aa.temp)
      bb <- rbind(bb, bb.temp)

      }
    
          aa <- aa[order(aa$filename, as.numeric(aa$time_stamp_sec)), ]
          bb <- bb[order(bb$filename, as.numeric(bb$time_stamp_sec)), ]
          out <- list(aa, bb)
          names(out) <-c ("frames", "whales")
          
    return(out)
      
     
}
      
    