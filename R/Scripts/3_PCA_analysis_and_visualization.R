#' ---
#' Title: PCA on full raw dataset and visualization
#' Author: Ana Eguiguren
#' Date: May 2026
#' Description: This script takes all frames collected during and outside
#' interactions (n = 258), and carries out a PCA




#0. LOAD DATA AND LIBRARIES -----
rm(list = ls())


library(factoextra)
library(scales)
library(dplyr)
library(scico)     # for colour palettes

source("R/Scripts/2_Prep_Data_For_Analysis.R")


# 1. Run PCA-----
pca.frames <- dat %>%
  dplyr::select(NumW, sd_head, Rank, File, PVA,Mean_NN, mean_body_angle, Maj_axis)%>%
  prcomp(scale = T)

#explore outputs

summary(pca.frames)

# 3. Scree Plot-----

tiff("Figures/PCA_formation_frame_scree_plot.tiff", 
     width = 7, height = 4.8, units = "in", res = 300)

fviz_eig(pca.frames, addlabels = T)

dev.off()


# 4. EXTRACT SCORES -----
frame.inds<-get_pca_ind(pca.frames)
frame.pca.scores <- frame.inds$coord
dim(frame.pca.scores)


frame.pca.scores <- bind_cols(dat, frame.pca.scores)


#edit behavioural states 
frame.pca.scores <- frame.pca.scores %>%
  mutate(intensity_sp = ifelse(is.na(intensity_sp), "not.interacting", intensity_sp))


# save for correlations

saveRDS(frame.pca.scores, 
        "PCA_Scores_all_frames.rds")

# 5. EXTRACT LOADINGS -----
frame.vars <- get_pca_var(pca.frames)
loadings_df <- as.data.frame(frame.vars$coord)
loadings_df$varname <- rownames(loadings_df)


write.csv(pca.frames$rotation, "Figures/PCA_loadings.csv")


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









#5. Biplots-----

#~~~~a. Overall biplot-----


lims <- c(-1,1)

tiff("Figures/PCA_all_frames_sociality.tiff", 
     width = 8, height = 4.8, units = "in", res = 300)

ggplot(
  frame.pca.scores,
  aes(x = Dim.1, y = Dim.2)
) +
  
  # background points (small, no outline)
  geom_point(
    data = subset(frame.pca.scores,
                  intensity_sp == "not.interacting"),
    aes(colour = soc),
    shape = 16,
    size = 2,
    alpha = 0.7
  ) +
  
  # foreground interacting points
  geom_point(
    data = subset(frame.pca.scores,
                  intensity_sp != "not.interacting"),
    aes(fill = soc,
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

