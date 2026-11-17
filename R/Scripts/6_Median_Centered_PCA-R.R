#6_median_Centered_PCA
# do median-centered PCA (full), and MANOVA

rm(list = ls())
library(factoextra)
library(scales)
library(dplyr)
library(lubridate)
library(scico)
set.seed(1234)

#source("R/Scripts/5_Median_Centered_Data_Prep.R")
dat_cent <- readRDS("Data/Processed_Data/median_centered_full_data_form.rds")

#inspect
head(dat_cent)


#1. Run PCA -----
pca.frames <- dat_cent %>%
  dplyr::select(NumW, sd_head, Rank, File, PVA,Mean_NN, mean_body_angle, Maj_axis)%>%
  prcomp(scale = T)


summary(pca.frames)

# 3. Scree Plot-----

tiff("Figures/formation_frame_scree_plot_cent.tiff", 
     width = 7, height = 4.8, units = "in", res = 300)

fviz_eig(pca.frames, addlabels = T)

dev.off()


# 4. Extract Scores -----
frame.inds<-get_pca_ind(pca.frames)
frame.pca.scores <- frame.inds$coord

frame.pca.scores <- bind_cols(dat_cent, frame.pca.scores)

#flip sign in PC2-R for consistency with PCA (raw)
frame.pca.scores <- frame.pca.scores %>%
  mutate(Dim.2 = Dim.2 * -1)

#process loadings: 
frame.vars <- get_pca_var(pca.frames)
loadings_df <- as.data.frame(frame.vars$coord)
loadings_df$varname <- rownames(loadings_df)


#flip sign in PC2-R for consistency with PCA (raw)
loadings_df <- loadings_df %>%
  mutate(Dim.2 = Dim.2 * -1)

arrow_scale <- 4  # Adjust as needed
loadings_df <- loadings_df %>%
  mutate(Dim.1 = Dim.1 * arrow_scale,
         Dim.2 = Dim.2 * arrow_scale)
#give nice names:
loadings_df<- loadings_df %>%
  mutate(varname = c("NumW", "SDHead", "Rank", "File","PVA", "MeanNN", 
                     "MeanBA", "MajAxis"))


pca_var_explained <- pca.frames$sdev^2 / sum(pca.frames$sdev^2)
percent_var <- round(100 * pca_var_explained, 1)  # Round to 1 decimal

# Step 2: Create axis labels with % variance
x_lab <- paste0("PC1 (", percent_var[1], "%)")
y_lab <- paste0("PC2 (", percent_var[2], "%)")


#save loadings for results
raw_loadings <- as.data.frame(pca.frames$rotation)

raw_loadings <- raw_loadings %>% mutate(PC2 = PC2 *-1)
write.csv(raw_loadings, "Figures/pca_loadings_all_cent.csv")




#5. Biplot-----
lims <- c(-1,1)
frame.pca.scores <- frame.pca.scores %>%
  mutate(intensity_sp = ifelse(
    is.na(intensity_sp), 
    "not.interacting", 
    intensity_sp
  ))


tiff("Figures/PCA-R_all_frames_sociality.tiff", 
     width = 8, height = 4.8, units = "in", res = 300)

ggplot(
  frame.pca.scores,
  aes(x = Dim.1, y = Dim.2)
) +
  
  # background points (small, no outline)
  geom_point(
    data = subset(frame.pca.scores,
                  intensity_sp == "not.interacting"),
    aes(colour = soc_raw),
    shape = 16,
    size = 2,
    alpha = 0.7
  ) +
  
  # foreground interacting points
  geom_point(
    data = subset(frame.pca.scores,
                  intensity_sp != "not.interacting"),
    aes(fill = soc_raw,
        shape = factor(intensity_sp)),
    size = 4.5,
    stroke = 0.8,
    colour = "black",
    alpha = 0.9
  ) +
  
  # axes
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black"
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "black"
  ) +
  
  # loading arrows
  geom_segment(
    data = loadings_df,
    aes(x = 0, y = 0,
        xend = Dim.1, yend = Dim.2),
    arrow = arrow(length = unit(0.25, "cm")),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.5
  ) +
  
  # variable labels
  geom_text(
    data = loadings_df,
    aes(x = Dim.1,
        y = Dim.2,
        label = varname),
    inherit.aes = FALSE,
    color = "black",
    vjust = -1,
    size = 4,
    hjust = 0.5
  ) +
  
  labs(
    x = x_lab,
    y = y_lab,
    shape = "Interaction intensity",
    fill = "Sociality"
  ) +
  
  # identical palettes for fill + colour
  scale_fill_scico(
    palette = "roma",
    limits  = c(lims[1], lims[2]),
    oob     = scales::squish,
    name    = "Sociality"
  ) +
  
  scale_colour_scico(
    palette = "roma",
    limits  = c(lims[1], lims[2]),
    oob     = scales::squish,
    guide = "none"
  ) +
  
  # foreground shapes only
  scale_shape_manual(
    values = c(
      "Aversive_Gm" = 21,
      "Aversive_Oo" = 24,
      "Neutral"     = 22
    ),
    name = "interaction state"
  ) +
  
  theme_classic()


dev.off()


#6. SAVE FOR GLMMS-----
saveRDS(frame.pca.scores, "Data/Processed_Data/PCA-R_median_centered_full_data_form.rds")


#**************************ARCHIVAL STUFF-----
#reorder factor levels 

risk_cols <- c("#721F81FF", "#CD4071FF",  "#FD9567FF")


ggplot(frame.pca.scores, 
             aes(x = Dim.1, y = Dim.2, 
                 colour = intensity_sp, 
                 shape = intensity_sp)) +
  # PCA points
  geom_point( size = 3, alpha = 0.7) +
  
  # Axes lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  
  # Arrows for loadings
  geom_segment(data = loadings_df,
               aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2),
               arrow = arrow(length = unit(0.25, "cm")),
               inherit.aes = FALSE, color = "black", linewidth = 0.5) +
  
  # Variable labels
  geom_text(data = loadings_df,
            aes(x = Dim.1, 
                y = Dim.2, 
                label = varname),
            inherit.aes = FALSE, 
            color = "black", 
            size = 3) +
  
  labs(x = x_lab, y = y_lab) +
  scale_color_manual(values = c("Aversive_Gm" = risk_cols[2],
                                "Neutral" = risk_cols[3], 
                                "Aversive_Oo" = risk_cols[1], 
                                "not.interacting" = "gray27"))+
  scale_shape_manual(values = c("Aversive_Gm" = 17,
                                "Neutral" = 17, 
                                "Aversive_Oo" = 17, 
                                "not.interacting" = 1))+
  
  theme_classic()


p.v

# 7. Plot Median centered-PCA values by event/aversity-----

#~~~~a. subset before/during/after frames-----
#format time
frame.varimax.scores <- frame.varimax.scores %>%
  mutate(frame_local_time = as.POSIXct(frame_local_time, format = "%Y-%m-%d %H:%M:%S") )

#get event timing

e_raw <- read.csv("Data/Processed_Data/Boris_Event_Times.csv", header = T)

#get start and end of each interaction (focal follow based)

e <- e_raw %>%
  group_by(event) %>%
  summarize(start_date_time = first(start_date_time), 
            end_date_time = last(end_date_time), 
            duration_focal = sum(duration)) %>%
  mutate(start_date_time = as.POSIXct(start_date_time, format = "%Y-%m-%d %H:%M:%S"), 
         end_date_time = as.POSIXct(end_date_time, format = "%Y-%m-%d %H:%M:%S") )

# join to event intensity and spp
# read in event descriptions:


e <- e %>%
  mutate(event.code = tolower(event))


d <- d %>%
  mutate(
    event.code = ifelse(
      grepl("_", event.code),
      event.code,
      gsub(" ", "_", event.code)
    )
  )

# join them

e <- left_join(e, d, by = 'event.code')




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


#subset observations during and close to events:
obs_near_events <- frame.varimax.scores %>%
  inner_join(
    event_windows,
    by = join_by(
      frame_local_time >= start_include,
      frame_local_time <= end_include
    )
  )




# filter out frames too close to event

obs_filtered <- obs_near_events %>%
  filter(
    !(
      frame_local_time >= start_exclude &
        frame_local_time <= end_exclude &
        (
          frame_local_time < start_date_time |
            frame_local_time > end_date_time
        )
    )
  )



# lable as before/during/after
obs_filtered <- obs_filtered %>%
  mutate(
    phase = case_when(
      frame_local_time >= start_date_time & frame_local_time <= end_date_time ~ "during",
      frame_local_time <  start_date_time ~ "before",
      frame_local_time >  end_date_time   ~ "after"
    )
  )




# create a sp. intensity column

obs_filtered <- obs_filtered %>%
  mutate(sp_intensity = paste(species.y, Intensity.y, sep = "_"))




# get centroids:

cent_obs <- obs_filtered %>%
  group_by(event.code.y, phase) %>%
  summarize(
    intensity_sp = first(intensity_sp.y), 
    first_frame = first(frame_local_time),
    last_frame = last(frame_local_time),
    vmax.1 = mean(vmax.1, na.rm = T),
    vmax.2 = mean(vmax.2, na.rm = T),
    NumW = mean(NumW, na.rm = T), 
    sd_head = mean(sd_head, na.rm = T), 
    Rank = mean(Rank, na.rm = T), 
    File = mean(File, na.rm = T),
    PVA = mean(PVA, na.rm = T), 
    Mean_NN = mean(Mean_NN, na.rm = T), 
    mean_body_angle = mean(mean_body_angle, na.rm = T),
    Maj_axi = mean(Maj_axis, na.rm = T)
  )%>%
  arrange(first_frame)


# save these for manovas:

obs_filtered <- obs_filtered %>%
  mutate(intensity_sp = intensity_sp.y)


# keep only events with at least 'during' state:
evs <- cent_obs %>%
  filter(phase == "during") %>% 
  dplyr::select(event.code.y, intensity_sp)%>%
  arrange(intensity_sp)


cent_obs <- cent_obs %>%
  filter(event.code.y %in% evs$event.code.y)

obs_filtered <- obs_filtered %>%
  filter(event.code.y %in% evs$event.code.y)



# keep only the ones with at least 2 states
obs_filtered <- obs_filtered %>%
  mutate(treat = ifelse(event_state =="not", "not", "during"), 
         intensity_sp = factor(intensity_sp, 
                               levels = c(
                                 "Neutral",
                                 "Aversive_Gm", 
                                 "Aversive_Oo"
                               )))




evs <- cent_obs %>%
  group_by(event.code.y)%>%
  filter(n()>1)%>%
  dplyr::select(event.code.y)%>%
  distinct()


obs_filtered <- obs_filtered %>%
  filter(event.code.y %in% evs$event.code.y)

cent_obs <- cent_obs %>%
  filter(event.code.y %in% evs$event.code.y)

levels(factor(cent_obs$intensity_sp))

# ~~~~b. plot all before/during/after frames with centroids-----

risk_cols <- c("#721F81FF", "#CD4071FF",  "#FD9567FF")



# extend axis limits for  labels:

pad <- 0.05  # 15% padding

xlim <- range(obs_filtered$vmax.1)
ylim <- range(obs_filtered$vmax.2)



# all (for legend:)



pA<- ggplot(cent_obs,
            aes(x = vmax.1, y = vmax.2,
                shape = phase, 
                colour = intensity_sp, 
                fill = intensity_sp))+
  
  
  #all points (before/during/after)
  geom_point(inherit.aes = F, 
             data = obs_filtered,
             aes(x = vmax.1, y = vmax.2, 
                  
                 shape = phase), 
             alpha = 0.5, 
             size = 2,
             colour = "grey")+
  #centroid points
  geom_point(size = 5, alpha = 0.8)+
  
  # Arrows for loadings
  geom_segment(data = v_loadings_df,
               aes(x = 0, y = 0, xend = vmax.1, yend = vmax.2),
               arrow = arrow(length = unit(0.25, "cm")),
               inherit.aes = FALSE, color = "black", linewidth = 0.5) +
  
  # Variable labels
  geom_text(data = v_loadings_df,
            aes(x = vmax.1, 
                y = vmax.2, 
                label = varname,
                vjust = vjust_lab,
                hjust = hjust_lab),
            inherit.aes = FALSE, 
            color = "black", 
            size = 4) +
  
  #reference lines:
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
  
  #colour setup
  scale_color_manual(values = c("Aversive_Gm" = risk_cols[2],
                                "Neutral" = risk_cols[3], 
                                "Aversive_Oo" = risk_cols[1]))+
  scale_fill_manual(values = c("Aversive_Gm" = risk_cols[2],
                               "Neutral" = risk_cols[3], 
                               "Aversive_Oo" = risk_cols[1]))+
  
  #shape setup
  scale_shape_manual(values = c("after" = 22, 
                                "before" = 21, 
                                "during" = 24))+
  
  #labels and style
  labs(x = x_lab, y = y_lab, title = "a)") +
  theme_classic()+
  coord_cartesian(
    xlim = xlim + c(-1, 1) * diff(xlim) * pad,
    ylim = ylim + c(-1, 1) * diff(ylim) * pad
  )+
  theme(legend.position = "null")



pA


ggsave("Figures/formation_frame_centroid_vmax_intensity_median_centered.png",
       pA, width = 5, height = 3.75)


#~~~~c. plot by sp_int-----


#~~~~~~~i. KW----
#make polygon
hull_KW <- obs_filtered %>%
  filter(
    event.code.y == "event_2",
    phase == "during"
  ) %>%
  slice(chull(vmax.1, vmax.2))


pKW<- ggplot(cent_obs %>% 
               filter(event.code.y == "event_2"), 
             aes(x = vmax.1, y = vmax.2,
                 shape = phase, 
                 colour = intensity_sp, 
                 fill = intensity_sp))+
  
  
  #centroid points
  geom_point(size = 4, alpha = 0.8)+
  
  #polygon
  geom_polygon(
    data = hull_KW,
    aes(vmax.1, vmax.2),
    fill = risk_cols[1],
    linewidth = 0.2,
    color = risk_cols[1],alpha = 0.2
  )+
  
  #all points (before/during/after)
  geom_point(inherit.aes = F, 
             data = obs_filtered %>%
               filter(event.code.y == "event_2"), 
             aes(x = vmax.1, y = vmax.2, 
                 colour  = intensity_sp, 
                 shape = phase), 
             alpha = 0.5
             ,size = 2)+
  
  #reference lines:
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
  
  #colour setup
  scale_color_manual(values = c("Aversive_Gm" = risk_cols[2],
                                "Neutral" = risk_cols[3], 
                                "Aversive_Oo" = risk_cols[1]))+
  scale_fill_manual(values = c("Aversive_Gm" = risk_cols[2],
                               "Neutral" = risk_cols[3], 
                               "Aversive_Oo" = risk_cols[1]))+
  
  #shape setup
  scale_shape_manual(values = c("after" = 22, 
                                "before" = 21, 
                                "during" = 24))+
  
  #labels and style
  labs(x = x_lab, y = y_lab, title = "b)") +
  theme_classic()+
  coord_cartesian(
    xlim = xlim + c(-1, 1) * diff(xlim) * pad,
    ylim = ylim + c(-1, 1) * diff(ylim) * pad
  )+
  theme(legend.position = "null")


pKW



ggsave("Figures/KW_formation_frame_centroid_vmax_intensity_med_cent.png",
       pKW, width = 4, height = 3)


#~~~~~~~~i. PW - Aversive-----

#make polygon
hull_PW <- obs_filtered %>%
  filter(
    intensity_sp == "Aversive_Gm" &
      phase == "during"
  ) %>%
  slice(chull(vmax.1, vmax.2))




pGmA<- ggplot(cent_obs %>% 
                filter(intensity_sp == "Aversive_Gm"), 
              aes(x = vmax.1, y = vmax.2,
                  shape = phase, 
                  colour = intensity_sp, 
                  fill = intensity_sp))+
  
  
  
  geom_point(size = 4, alpha = 0.8)+
  
  #polygon
  geom_polygon(
    data = hull_PW,
    aes(vmax.1, vmax.2),
    fill = risk_cols[2],
    linewidth = 0.2,
    color = risk_cols[2],alpha = 0.2
  )+
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray", 
             linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray", 
             linewidth = 0.5) +
  
  
  geom_point(inherit.aes = F, 
             data = obs_filtered %>%
               filter(intensity_sp == "Aversive_Gm"), 
             aes(x = vmax.1, y = vmax.2, 
                 colour  = intensity_sp, 
                 shape = phase), 
             alpha = 0.5, 
             size = 2)+
  
  scale_color_manual(values = c("Aversive_Gm" = risk_cols[2],
                                "Neutral" = risk_cols[3], 
                                "Aversive_Oo" = risk_cols[1]))+
  scale_fill_manual(values = c("Aversive_Gm" = risk_cols[2],
                               "Neutral" = risk_cols[3], 
                               "Aversive_Oo" = risk_cols[1]))+
  
  
  scale_shape_manual(values = c("after" = 22, 
                                "before" = 21, 
                                "during" = 24))+
  
  labs(x = x_lab, y = y_lab, title = "c)") +
  
  theme_classic()+
  
  coord_cartesian(
    xlim = xlim + c(-1, 1) * diff(xlim)*pad ,
    ylim = ylim + c(-1, 1) * diff(ylim)*pad 
  )+
  theme(legend.position = "null")


pGmA


ggsave("Figures/GmA_formation_frame_centroid_vmax_intensity_med_cent.png",
       pGmA, width = 4, height = 3)



#~~~~~~~~i. Neutral-----
#make polygon
hull_N <- obs_filtered %>%
  filter(
    intensity_sp == "Neutral" &
      phase == "during"
  ) %>%
  slice(chull(vmax.1, vmax.2))


pN<- 
  ggplot(cent_obs %>% 
           filter(intensity_sp == "Neutral" ), 
         aes(x = vmax.1, y = vmax.2,
             shape = phase, 
             colour = intensity_sp, 
             fill = intensity_sp))+
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray",
             linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray", 
             linewidth = 0.5) +
  
  
  geom_point(size = 4, alpha = 0.8)+
  
  #polygon
  geom_polygon(
    data = hull_N,
    aes(vmax.1, vmax.2),
    fill = risk_cols[3],
    linewidth = 0.2,
    color = risk_cols[3],alpha = 0.2
  )+
  
  geom_point(inherit.aes = F, 
             data = obs_filtered %>%
               filter(intensity_sp == "Neutral"), 
             aes(x = vmax.1, y = vmax.2, 
                 colour  = intensity_sp, 
                 shape = phase), 
             alpha = 0.4, 
             size = 2)+
  
  
 
  scale_color_manual(values = c("Aversive_Gm" = risk_cols[2],
                                "Neutral" = risk_cols[3], 
                                "Aversive_Oo" = risk_cols[1]))+
  scale_fill_manual(values = c("Aversive_Gm" = risk_cols[2],
                               "Neutral" = risk_cols[3], 
                               "Aversive_Oo" = risk_cols[1]))+
  
  scale_shape_manual(values = c("after" = 22, 
                                "before" = 21, 
                                "during" = 24))+
  labs(x = x_lab, y = y_lab, title = "d)") +
  theme_classic()+
  coord_cartesian(
    xlim = xlim + c(-1, 1) * diff(xlim) * pad,
    ylim = ylim + c(-1, 1) * diff(ylim) * pad
  )+
  theme(legend.position = "null")


pN


ggsave("Figures/pN_formation_frame_centroid_vmax_intensity_med_cent.png",
       pN, width = 4, height = 3)



#~~~~~~ save plot -----
pa <- (pA + theme(legend.position = "Null"))

p <- pKW / pGmA / pN


pa + p






ggsave("Figures/PCA_vmax_centroids_during_not_med_cent_all.png",
       p, width =3.5, height = 8)


# 8. PCs Violin Plots -----

# keep only the ones with at least 2 states
obs_filtered <- obs_filtered %>%
  mutate(treat = ifelse(event_state =="not", "not", "during"), 
         intensity_sp = factor(intensity_sp, 
                               levels = c(
                                 "Neutral",
                                 "Aversive_Gm", 
                                 "Aversive_Oo"
                               )))



head(data.frame(cent_obs))

evs <- cent_obs %>%
  group_by(event.code.y)%>%
  filter(n()>1)%>%
  dplyr::select(event.code.y)%>%
  distinct()


obs_filtered <- obs_filtered %>%
  filter(event.code.y %in% evs$event.code.y)



pc1<- ggplot(
  data = obs_filtered, 
  aes(
    x = intensity_sp,
    y = vmax.1, 
    colour = interaction(treat, intensity_sp), 
    fill = interaction(treat, intensity_sp)
  )
) +
  
  geom_hline(yintercept = 0, 
             linetype = "dashed")+
  
  geom_boxplot(
    position = position_dodge(width = 1),
    outlier.shape = NA, 
    alpha = 0.5
  ) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.5,
      dodge.width  = 1
    ),
    alpha = 0.7, 
    size = 2
  ) +
  facet_wrap(~intensity_sp,
             scales = "free_x")+
  scale_colour_manual(
    values = c(
      "not.Neutral"        = "grey60",
      "not.Aversive_Gm"    = "grey60",
      "not.Aversive_Oo"    = "grey60",
      "during.Neutral"     = risk_cols[3],
      "during.Aversive_Gm" = risk_cols[2],
      "during.Aversive_Oo" = risk_cols[1]
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "not.Neutral"        = "grey60",
      "not.Aversive_Gm"    = "grey60",
      "not.Aversive_Oo"    = "grey60",
      "during.Neutral"     = risk_cols[3],
      "during.Aversive_Gm" = risk_cols[2],
      "during.Aversive_Oo" = risk_cols[1]
    )
  ) +
  
  labs(x = "", 
       y = "PC1 Score",
       title = "a)") +
  
  theme_classic()+
  theme(legend.position = "")

pc1 <-pc1 + theme(strip.text = element_blank())





pc2<- ggplot(
  data = obs_filtered, 
  aes(
    x = intensity_sp,
    y = vmax.2, 
    colour = interaction(treat, intensity_sp), 
    fill = interaction(treat, intensity_sp)
  )
) +
  
  geom_hline(yintercept = 0, 
             linetype = "dashed")+
  
  geom_boxplot(
    position = position_dodge(width = 1),
    outlier.shape = NA, 
    alpha = 0.5
  ) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.5,
      dodge.width  = 1
    ),
    alpha = 0.7, 
    size = 2
  ) +
  facet_wrap(~intensity_sp,
             scales = "free_x")+
  scale_colour_manual(
    values = c(
      "not.Neutral"        = "grey60",
      "not.Aversive_Gm"    = "grey60",
      "not.Aversive_Oo"    = "grey60",
      "during.Neutral"     = risk_cols[3],
      "during.Aversive_Gm" = risk_cols[2],
      "during.Aversive_Oo" = risk_cols[1]
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "not.Neutral"        = "grey60",
      "not.Aversive_Gm"    = "grey60",
      "not.Aversive_Oo"    = "grey60",
      "during.Neutral"     = risk_cols[3],
      "during.Aversive_Gm" = risk_cols[2],
      "during.Aversive_Oo" = risk_cols[1]
    )
  ) +
  
  labs(x = "", 
       y = "PC2 Score", 
       title = "b)") +
  
  theme_classic()+
  theme(legend.position = "none")

pc2 <- pc2+ theme(strip.text = element_blank())




bps <- pc1 + pc2
bps

# save:
ggsave("Figures/Boxplots_PCs_during_not_median_centered.png",
       bps, width = 8, height = 3.5)


# 9. PERMANOVA-----

# only big rule is homogeneity of variances 
library(heplots)
#make a grouping variable: 


#~~~~a. test assumptions------
# homogeneity of variances (normality is not important)
# should scale 
var_dep <- data.frame(Dim1 = scale(obs_filtered$vmax.1), 
                      Dim2 = scale(obs_filtered$vmax.2))


hist(var_dep$Dim1, breaks = 20)
hist(var_dep$Dim2, breaks = 20)



obs_filtered <- obs_filtered %>%
  mutate(grouping = paste(treat, intensity_sp, sep = "_"))


boxM(var_dep, factor(obs_filtered$grouping))

#absolutely broken!

# evennes of distributions?
library(vegan)

perm_dist <- vegdist(var_dep, method = 'euclidean')

dispersion <- betadisper(perm_dist, group = obs_filtered$grouping)


# site scores (PCoA coordinates)
sites <- as.data.frame(scores(dispersion, display = "sites"))
sites$group <- dispersion$group

# centroid scores
centroids <- as.data.frame(scores(dispersion, display = "centroids"))
centroids$group <- rownames(centroids)


#convex hulls
hulls <- sites %>%
  group_by(group) %>%
  filter(n() >= 3) %>%
  slice(chull(PCoA1, PCoA2)) %>%
  reframe(
    PCoA1 = c(PCoA1, first(PCoA1)),
    PCoA2 = c(PCoA2, first(PCoA2)),
    group = c(group, first(group))
  )


ggplot(sites, aes(PCoA1, PCoA2, colour = group)) +
  geom_polygon(data = hulls,
               aes(fill = group),
               alpha = 0.1,
               colour = NA) +
  geom_path(data = hulls,
            aes(group = group),
            linewidth = 0.8) +
  #geom_point(size = 2) +
  geom_point(data = centroids,
             aes(PCoA1, PCoA2, colour = group),
             shape = 4, size = 5, stroke = 1.2) +
  coord_equal() +
  theme_classic() +
  labs(colour = "Group", fill = "Group")


#unequal dispersion!
anova(dispersion) #super significant

#~~~~b. run permanova----
pman_pc<- adonis2(perm_dist ~ intensity_sp*treat, data = obs_filtered, 
                  permutations = 999, by = "margin", 
                  strata = obs_filtered$date)



print(pman_pc)#
#interaction term is not significant

pman_pc2<- adonis2(perm_dist ~ intensity_sp*treat, data = obs_filtered, 
                  permutations = 999, by = "term")

print(pman_pc2)#

obs_filtered$dist <- dispersion$distances


ggplot(obs_filtered,
       aes(treat, dist, fill = treat)) +
  geom_boxplot() +
  facet_wrap(~ intensity_sp) +
  theme_bw(base_size = 10)








#10. Centroid-based PERMANOVA-----

mean_obs<- obs_filtered %>%
  group_by(event.code.y, treat) %>%
  summarize(n = n(), 
            sp_intensity = first(intensity_sp.y),
            treat = first(treat), 
            Dim1 = mean(vmax.1, na.rm = T), 
            Dim2 = mean(vmax.2, na.rm = T), 
            date = first(date)
            )



#~~~~a. check assumptions----

# homogeneity of variances (normality is not important)
# should scale 
var_dep_m <- data.frame(Dim1 = scale(mean_obs$Dim1), 
                      Dim2 = scale(mean_obs$Dim2))


hist(var_dep_m$Dim1, breaks = 20)
hist(var_dep_m$Dim2, breaks = 20)



mean_obs <- mean_obs %>%
  mutate(grouping = paste(treat, sp_intensity, sep = "_"))


boxM(var_dep_m, factor(mean_obs$grouping))

#too few data points
# evennes of distributions?

perm_dist_m <- vegdist(var_dep_m, method = 'euclidean')

dispersion <- betadisper(perm_dist_m, group = mean_obs$grouping)


# site scores (PCoA coordinates)
sites <- as.data.frame(scores(dispersion, display = "sites"))
sites$group <- dispersion$group

# centroid scores
centroids <- as.data.frame(scores(dispersion, display = "centroids"))
centroids$group <- rownames(centroids)


#convex hulls
hulls <- sites %>%
  group_by(group) %>%
  filter(n() >= 3) %>%
  slice(chull(PCoA1, PCoA2)) %>%
  reframe(
    PCoA1 = c(PCoA1, first(PCoA1)),
    PCoA2 = c(PCoA2, first(PCoA2)),
    group = c(group, first(group))
  )


ggplot(sites, aes(PCoA1, PCoA2, colour = group)) +
  geom_polygon(data = hulls,
               aes(fill = group),
               alpha = 0.1,
               colour = NA) +
  geom_path(data = hulls,
            aes(group = group),
            linewidth = 0.8) +
  #geom_point(size = 2) +
  geom_point(data = centroids,
             aes(PCoA1, PCoA2, colour = group),
             shape = 4, size = 5, stroke = 1.2) +
  coord_equal() +
  theme_classic() +
  labs(colour = "Group", fill = "Group")


#unequal dispersion!
anova(dispersion) #super significant

#~~~~b. run permanova----
pman_pc_m<- adonis2(perm_dist_m ~ sp_intensity*treat, data = mean_obs, 
                  permutations = 999, by = "term")



print(pman_pc_m)#
#interaction term is not significant

print(pman_pc2)#

mean_obs$dist <- dispersion$distances


ggplot(mean_obs,
       aes(treat, dist, fill = treat)) +
  geom_boxplot() +
  facet_wrap(~ sp_intensity) +
  theme_bw(base_size = 10)




#~~~~~c. Centroid based MANOVA----
var_dep_t <- t(var_dep_m)

mshapiro.test(var_dep_t) # non normal (tiny dataset)


MAN_cent_var <- manova(as.matrix(var_dep_m) ~ sp_intensity*treat, 
                       data = mean_obs)
summary(MAN_cent_var)


# GLMM for each PC-----

# ~~~a. PC1----

summary(factor(obs_filtered$intensity_sp))
#because it becomes singular, try out lumping kws and gms?
obs_filtered <- obs_filtered %>%
  mutate(intensity = ifelse(intensity_sp == "Neutral", "Neutral", "Aversive"))


summary(factor(obs_filtered$intensity))

library(nlme)
library("AICcmodavg")




m1.lmer <- lme(vmax.1 ~ treat * intensity, 
               random = ~1 |date, 
                 data = obs_filtered, 
               method = "ML")


AICc(m1.lmer, second.ord = T)


m1.lm <- gls(vmax.1 ~ treat * intensity, 
             method = "ML", data = obs_filtered)

AICc(m1.lm, second.ord = T)




acf(resid(m1.lm)) # there is dependence here!
acf(resid(m1.lmer))

#save acf
# compute ACF objects without plotting
acf_lm   <- acf(resid(m1.lm),   plot = FALSE)
acf_lmer <- acf(resid(m1.lmer), plot = FALSE)

# find common y-axis limits
ylim_common <- range(
  acf_lm$acf,
  acf_lmer$acf
)

# save the plot
png("Figures/Model_Diagnostics/acf_comparison.png",
    width = 1600,   # double the pixels
    height = 800,
    res = 300)

par(mfrow = c(1, 2),
    mar = c(4, 2, 1, 1))  # two panels side by side

plot(acf_lm,
     ylim = ylim_common,
     main = "a)")

plot(acf_lmer,
     ylim = ylim_common,
     main = "b)",
     yaxt = "n",   # no y-axis ticks
     ylab = "") 

dev.off()


anova(m1.lmer, m1.lm)# random structure does help 


plot(m1.lm)
plot(factor(obs_filtered$treat), resid(m1.lm))
plot(factor(obs_filtered$intensity), resid(m1.lm))
plot(interaction(obs_filtered$intensity,
                   obs_filtered$treat), resid(m1.lm))#looks ok!



summary(m1.lmer)

# explore keeping variables


m1.lm.a <- lme(vmax.1 ~ treat + intensity, 
               random = ~1 |date, 
               data = obs_filtered, 
               method = "ML")

AICc(m1.lm.a)

anova(m1.lmer, m1.lm.a) # removing interaction is good



m1.lm.b <- lme(vmax.1 ~ treat , 
               random = ~1 |date, 
               data = obs_filtered, 
               method = "ML")

m1.lm.c <- lme(vmax.1 ~  intensity, 
               random = ~1 |date, 
               data = obs_filtered, 
               method = "ML")


AICc(m1.lm.b)
AICc(m1.lm.c)


anova(m1.lm.a, 
      m1.lm.b, 
      m1.lm.c)


# intensity alone is best


m1.lm.0 <- lme(vmax.1 ~  1, 
               random = ~1 |date, 
               data = obs_filtered, 
               method = "ML",
               weights = varIdent(form = ~1 | intensity))


AICc(m1.lm.0)
plot(m1.lm.0)

anova(m1.lm.0, m1.lm.c) # empty is best!

summary(m1.lm.0)




plot(m1.lm.c)


durbinWatsonTest(m1.lm.0)



# ~~~a. PC2----


m2.lmer <- lme(vmax.2 ~ treat * intensity, 
                random = ~1 | date, 
                data = obs_filtered, 
                method = "ML")




summary(m2.lmer)


AICc(m2.lmer)

m2.lm <- gls(vmax.2 ~ treat * intensity, 
             method = "ML", data = obs_filtered)

AICc(m2.lm)



anova(m2.lmer, m2.lm)
acf(resid(m2.lm)) # tiny bit


acf(resid(m2.lmer))



# compute ACF objects without plotting
acf_lm_2   <- acf(resid(m2.lm),   plot = FALSE)
acf_lmer_2 <- acf(resid(m2.lmer), plot = FALSE)

# find common y-axis limits
ylim_common <- range(
  acf_lm_2$acf,
  acf_lm_2$acf
)

# save the plot
png("Figures/Model_Diagnostics/acf_comparison_PC2.png",
    width = 1600,   # double the pixels
    height = 800,
    res = 300)

par(mfrow = c(1, 2),
    mar = c(4, 2, 1, 1))  # two panels side by side

plot(acf_lm_2,
     ylim = ylim_common,
     main = "a)")

plot(acf_lmer_2,
     ylim = ylim_common,
     main = "b)",
     yaxt = "n",   # no y-axis ticks
     ylab = "") 

dev.off()




plot(m2.lm)
plot(factor(obs_filtered$treat), resid(m2.lm))
plot(factor(obs_filtered$intensity), resid(m2.lm))
plot(interaction(obs_filtered$intensity,
                 obs_filtered$treat), resid(m2.lm))


m2.lm.a <- gls(vmax.2 ~ treat + intensity, 
               method = "ML", data = obs_filtered)


AICc(m2.lm.a)



m2.lm.b <- gls(vmax.2 ~ treat, 
               method = "ML", data = obs_filtered)


AICc(m2.lm.b)


m2.lm.c <- gls(vmax.2 ~ intensity, 
               method = "ML", data = obs_filtered)


AICc(m2.lm.c)



m2.lm.0 <- gls(vmax.2 ~ 1, 
               method = "ML", data = obs_filtered)


AICc(m2.lm.0)




anova(m2.lm, m2.lm.a)

summary(m2.lm)
anova(m2.lm)
summary(m2.lm.a)


m2.lm$coefficients

#Value Std.Error   t-value p-value
#(Intercept)               -1.292192 0.3558479 -3.631305  0.0004
#treatnot                   1.608087 0.3907009  4.115905  0.0001
#intensityNeutral           1.903837 0.8716459  2.184187  0.0309
#treatnot:intensityNeutral -2.022322 0.9174155 -2.204369  0.0294

#very similar std error, p-values change a bit (random effects model is a bit 
#more conservative)


m2.lmer$coefficients # these are very similar, but the random effects model 
# drops the interaction term. 



library(emmeans)



emmeans(m2.lm, ~ treat | intensity)


summary(m2.lm)
anova(m2.lm)
# 
m2.lmer


m2.lmer.a <- lme(vmax.2 ~ treat + intensity, 
               random = ~1 | date, 
               data = obs_filtered, 
               method = "ML")


anova(m2.lmer.a, m2.lmer) # BIC same, AIC diff (?)some support! 



contr_m2 <- contrast(
  emm_m2,
  method = "pairwise",
  by = "intensity_sp",
  adjust = "tukey"
)
contr_m2 <- as.data.frame(contr_m2)





contr_m2$signif <- cut(
  contr_m2$p.value,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "ns")
)







# ~~~~a. plot----
library(emmeans)
emm_m1 <- emmeans(m1.lmer, ~ treat * intensity_sp)
emm_m2 <- emmeans(m2.lmer, ~ treat * intensity_sp)


emm_m1_df <- as.data.frame(emm_m1)
emm_m2_df <- as.data.frame(emm_m2)


# make significance things:
sig_df <- contr_m1 %>%
  mutate(
    y.position = max(emm_m1_df$upper.CL) * 1.05
  )

sig_df2 <- contr_m2 %>%
  mutate(
    y.position = max(emm_m2_df$upper.CL) * 1.05
  )


ggplot(emm_m1_df,
       aes(x = intensity_sp, y = emmean, color = treat)) +
  geom_point(
    position = position_dodge(width = 0.4),
    size = 3
  ) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.2,
    position = position_dodge(width = 0.4)
  ) +
  
  geom_text(
    data = sig_df,
    aes(x = intensity_sp, y = y.position, label = signif),
    inherit.aes = FALSE,
    size = 5
  ) +
  
  labs(
    y = "Estimated PC1",
    x = "Intensity"
  ) +
  theme_classic()





ggplot(emm_m2_df,
       aes(x = intensity_sp, y = emmean, color = treat)) +
  geom_point(
    position = position_dodge(width = 0.4),
    size = 3
  ) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.2,
    position = position_dodge(width = 0.4)
  ) +
  
  geom_text(
    data = sig_df2,
    aes(x = intensity_sp, y = y.position, label = signif),
    inherit.aes = FALSE,
    size = 5
  ) +
  labs(
    y = "Estimated PC2",
    x = "Intensity"
  ) +
  theme_classic()


# still significant with bonferroni correction
# I sincerely trust this more! because it accounts for autocorrelation, 
# only does two tests, and I'm happy about it! :D bless this soul!!!



