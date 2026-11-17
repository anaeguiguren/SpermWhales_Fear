#Boxplots

#' 
#' ---
#' Title: Boxplots of PCA-R outcomes from cluster formation metrics 
#' Author: Ana Eguiguren
#' Date: May 2026
#' Description: Makes boxplots from subset of frames observed 
#' during/before/after events. 
#' ---



#0. LOAD LIBRARIES-----
rm(list = ls())

#list of files
packages <- c("dplyr", "tidyr", 
              "lubridate", "ggplot2")

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
library(patchwork)
library(ggplot2)


#1. LOAD DATA ----
#source("R/Scripts/7_Subset_PCA-R_for_GLMMs.R")

#obs_filtered <- readRDS("Data/Processed_Data/obs_filtered_med_cent.rds")

obs_filtered <- readRDS("Data/Processed_Data/obs_unfiltered_med_cent.rds")



obs_filtered$treat <- obs_filtered$interaction_window



# 2. CLEAN UP FOR PLOTTING-----

# Reorganize factor levels
# For intensity_sp: Neutral, Aversive_Gm, Aversive_Kw
obs_filtered$intensity_sp <- factor(obs_filtered$intensity_sp,
                                    levels = c("Neutral", "Aversive_Gm", "Aversive_Oo"))


#for behav_state

obs_filtered$behav_state_int <- factor(obs_filtered$behav_state_int,
                                    levels = c("Neutral", "Aversive_Gm", "Aversive_Oo"))


# For treat: 'more than 2 hours apart', 'up to two hours apart', 'during'
obs_filtered$treat <- factor(obs_filtered$treat,
                             levels = c("more than 2 hours", 
                                        "within 2hrs", 
                                        "during"))


#short lables
obs_filtered$treat_short <- dplyr::recode(
  obs_filtered$treat, 
  "more than 2 hours" = ">2hrs",
  "within 2hrs" = "<2hrs",
  "during" = "during"
)


obs_filtered$event_short <- factor(obs_filtered$event.id)


# prepare colours and labels-----
risk_cols <- c("#721F81FF", "#CD4071FF",  "#FD9567FF")

panel_labels <- c(
  "Neutral" = "Neutral",
  "Aversive_Gm" = "Aversive (GM)",
  "Aversive_Oo" = "Aversive (KW)"
)


# point size 

p.size = 2
# 3. Make plots ------
#~3a. PC1-R----
pc1 <- ggplot(
  data = obs_filtered, 
  aes(
    interaction(treat_short, event_short, sep = "\n"),    
    y = Dim.1, 
    colour = interaction(treat_short, intensity_sp), 
    fill = interaction(treat_short, intensity_sp)
  )
) +
  
  geom_hline(yintercept = 0, 
             linetype = "dashed") +
  
  geom_boxplot(
    outlier.shape = NA, 
    width = 0.6,
    alpha = 0.5
  ) +
  geom_jitter(
    width = 0.2,
    alpha = 0.5, 
    size = p.size
  ) +
  
  scale_colour_manual(
    values = c(
      # more than 2 hours apart - all grey
      "more than 2 hours.Neutral"     = "grey40",
      "more than 2 hours.Aversive_Gm" = "grey40",
      "more than 2 hours.Aversive_Oo" = "grey40",
      # up to two hours apart - all grey
      "within 2hrs.Neutral"       = "grey40",
      "within 2hrs.Aversive_Gm"   = "grey40",
      "within 2hrs.Aversive_Oo"   = "grey40",
      # during - different color for each intensity
      "during.Neutral"                      = risk_cols[3],
      "during.Aversive_Gm"                  = risk_cols[2],
      "during.Aversive_Oo"                  = risk_cols[1]
    )
  ) +
  
  scale_fill_manual(
    values = c(
      # more than 2 hours apart - all grey
      "more than 2 hours.Neutral"     = NA,
      "more than 2 hours.Aversive_Gm" = NA,
      "more than 2 hours.Aversive_Oo" = NA,
      # up to two hours apart - all grey
      "within 2hrs.Neutral"       = NA,
      "within 2hrs.Aversive_Gm"   = NA,
      "within 2hrs.Aversive_Oo"   = NA,
      # during - different color for each intensity
      "during.Neutral"                      = risk_cols[3],
      "during.Aversive_Gm"                  = risk_cols[2],
      "during.Aversive_Oo"                  = risk_cols[1]
    )
  ) +
  
  facet_wrap(~ behav_state_int, 
             nrow = 1,
             labeller = as_labeller(panel_labels),
             scales = "free_x", 
             space = "free_x") +
  
  labs(x = "", 
       y = "PC1-R Score",
       title = "a)") +
  
  theme_classic() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey90", colour = "black"),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
  )

#see
pc1

#~3b. PC2-R------
pc2 <- ggplot(
  data = obs_filtered, 
  aes(
    interaction(treat_short, event_short, sep = "\n"),    
    y = Dim.2, 
    colour = interaction(treat_short, intensity_sp), 
    fill = interaction(treat_short, intensity_sp)
  )
) +
  
  geom_hline(yintercept = 0, 
             linetype = "dashed") +
  
  geom_boxplot(
    outlier.shape = NA, 
    width = 0.6,
    alpha = 0.5
  ) +
  geom_jitter(
    width = 0.2,
    alpha = 0.5, 
    size = p.size
  ) +
  
  scale_colour_manual(
    values = c(
      # more than 2 hours apart - all grey
      "more than 2 hours.Neutral"     = "grey40",
      "more than 2 hours.Aversive_Gm" = "grey40",
      "more than 2 hours.Aversive_Oo" = "grey40",
      # up to two hours apart - all grey
      "within 2hrs.Neutral"       = "grey40",
      "within 2hrs.Aversive_Gm"   = "grey40",
      "within 2hrs.Aversive_Oo"   = "grey40",
      # during - different color for each intensity
      "during.Neutral"                      = risk_cols[3],
      "during.Aversive_Gm"                  = risk_cols[2],
      "during.Aversive_Oo"                  = risk_cols[1]
    )
  ) +
  
  scale_fill_manual(
    values = c(
      # more than 2 hours apart - all grey
      "more than 2 hours.Neutral"     = NA,
      "more than 2 hours.Aversive_Gm" = NA,
      "more than 2 hours.Aversive_Oo" = NA,
      # up to two hours apart - all grey
      "within 2hrs.Neutral"       = NA,
      "within 2hrs.Aversive_Gm"   = NA,
      "within 2hrs.Aversive_Oo"   = NA,
      # during - different color for each intensity
      "during.Neutral"                      = risk_cols[3],
      "during.Aversive_Gm"                  = risk_cols[2],
      "during.Aversive_Oo"                  = risk_cols[1]
    )
  ) +
  
  facet_wrap(~ behav_state_int, 
             nrow = 1,
             labeller = as_labeller(panel_labels),
             scales = "free_x", 
             space = "free_x") +
  
  labs(x = "", 
       y = "PC2-R Score",
       title = "b)") +
  
  theme_classic() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey90", colour = "black"),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
  )

#see

pc2

#~3c.make joint plot----
pc.boxplots <- pc1 / pc2
pc.boxplots
# save:
ggsave("Figures/Boxplots_PCs_during_not_by_event.png",
       pc.boxplots, width = 9, height = 6.5, 
       dpi = 300)
