#' ---
#' Title: Correlation between variables & sociality
#' Author: Ana Eguiguren
#' Date: May 2026
#' Description: This script uses the average formation metrics for each 
#' flight sessionto explore the relationships between formation metrics 
#' and sociality




#0. LOAD DATA AND LIBRARIES -----
rm(list = ls())
library(tidyr)
library(dplyr)
library(corrplot)



dat <- readRDS("PCA_Scores_all_frames.rds")


#1. AVERAGE PER FLIGHT SESSION FOR INDEPENDENCE -----

  
summary_table <- dat %>%
  group_by(flight_session_ID) %>%
  summarise(
    across(
      c(NumW, sd_head, Rank, File, PVA, Mean_NN, mean_body_angle, Maj_axis, 
        Dim.1, Dim.2, soc),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )


#2. PEARSON CORRELATIONS -----

#~2a. compute correlations----

# keep only the mean columns
mean_data <- summary_table %>%
  select(ends_with("_mean"))


# correlation matrix
cor_matrix <- cor(
  mean_data,
  use = "pairwise.complete.obs",
  method = "pearson"
)
rownames(cor_matrix) <- c("NumW", "SDHead", "Rank",
                          "File", "PVA", "MeanNN",
                          "MeanBA", "MajAxis", "PC1", 
                          "PC2", "Soc"
                          )

colnames(cor_matrix) <- rownames(cor_matrix)


#~2b. visualize correlations----

#method = c("circle", "square", "ellipse", "number", "shade", "color", "pie"),


corrplot(
  cor_matrix,
  method = "circle",
  type = "upper",
  order = "hclust",
  addCoef.col = "black",   # correlation coefficient color
  number.cex = 0.7,        # smaller correlation coefficient size
  number.font = 2,         # plain text (not bold)
  tl.col = "black",
  tl.cex = 0.8,            # label size
  tl.font = 1,             # plain text labels (not bold)
  tl.srt = 45,
  diag = FALSE
)

#3. SAVE CORRELATION MATRIX -----

#keep only lower half:

cor_matrix[upper.tri(cor_matrix)] <- NA

#~2b. SAVE-----


write.csv(cor_matrix, "Figures/Correlation_matrix_forms_soc.csv")
