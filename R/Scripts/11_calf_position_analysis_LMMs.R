
#' ------------------------------------------------------------------------------
#' Title: GLMM on calf position analysis-----
#' Author: Ana Eguiguren
#' Date: May 2026
#' Description:
#' ------------------------------------------------------------------------------

# Clear workspace

rm(list = ls())
#source("R/Scripts/10_calf_position_analysis_data_prep.R")
set.seed(12345)


library(dplyr)
library(patchwork)
library(lubridate)
library(ggplot2)
library(tidyr)
library(emmeans) #to extract effect sizes
library(MuMIn) # for model Selection
library(nlme) 


#1. LOAD DATA ----

#data by block
dat.block <- readRDS(
  "Data/Processed_Data/calf_position_block_summaries.RData"
)


#data by flight session
dat<- readRDS(
  "Data/Processed_Data/calf_position_individual_data.RData"
)

# 2. PREPARE DATA----

# summarize by cluster
dat.clust <- dat %>%
  group_by(cluster_ID, age_class_coarse)%>%
  summarise(
    date = first(date), 
    frame_ID = first(frame_ID), 
    flight_session_ID = first(flight_session_ID),
    phase = first(phase),
    intensity_sp = first(intensity_sp),
    NND = mean(NND, na.rm = T),
    filep_f_std = mean(filep_f_std, na.rm = T),
    rankP_min_std = mean(rankP_min_std, na.rm = T),
    dist_std = mean(dist_std, na.rm = T),
    prop_outer_25 = sum(outer_25)/ n(), 
    
    
    n_total = n(),
    n_outer = sum(outer_25, na.rm = TRUE),
    n_inner = n_total - n_outer
  )







# make a good state factor (i.e., during interaction vs. not)

dat.clust <- dat.clust %>%
  mutate(state = 
           factor(ifelse(is.na(intensity_sp), 
                         "Not",
                         ifelse(intensity_sp == "Neutral", 
                                "Not", 
                                "Aversive"))))

#~2b. general summaries----
dat.clust %>%
  filter(age_class_coarse == "Adult")%>%
  group_by(state)%>%
  summarize(n = n())%>%
  ungroup()

dat.clust$state <- factor(dat.clust$state)



# relevel (during interactions is the comparison level)
dat.clust$state <- relevel(dat.clust$state, ref = "Not")



# 3. MODEL SELECTION-----
# ~3a. NND -----
#~~~~ i.Select correct random structure (REML)----

# Null model (no random effects)
m1.ls <- gls(log(NND) ~ 1 + age_class_coarse * state,
             data = dat.clust,
             method = "REML")

summary(m1.ls)

# Model with random intercept for date
m1.ran <- lme(log(NND) ~ 1 + age_class_coarse * state,
              random = ~ 1 | flight_session_ID,
              data = dat.clust,
              method = "REML")


# Model with both random intercept and variance structure
m1.var.ran <- lme(log(NND) ~ 1 + age_class_coarse * state,
                  data = dat.clust,
                  random = ~ 1 | flight_session_ID,
                  weights = varIdent(form = ~ 1 | age_class_coarse),
                  method = "REML")



# save
m1.rand.tab <- tibble(
  Model = c(
    "No random structure",
    "Random intercept (flight session)",
    "Random intercept + Variance structure"
  ),
  
  AICc = c(
    AICc(m1.ls),
    AICc(m1.ran),
    AICc(m1.var.ran)
  )
) %>%
  
  mutate(
    deltaAICc = AICc - min(AICc)
  ) %>%
  
  arrange(AICc)

m1.rand.tab

#round to 2 decimal points
m1.rand.tab <- m1.rand.tab %>%
  mutate(
    across(c(AICc, deltaAICc),
           round, 2)
  )

write.csv(m1.rand.tab, "Figures/CalfPos_Model_Selection_Rand_Structure_MeanNND.csv")


#~~~~ ii.Select Fixed effects (ML)----
m1.final<- update(m1.ran, method = "ML")

# Run dredge and keep original object
m1.dredge_results <- dredge(m1.final)

# Convert to dataframe
m1.cand_tab_clean <- as.data.frame(m1.dredge_results)

# Add formulas
m1.cand_tab_clean$formula <- sapply(
  get.models(m1.dredge_results, subset = TRUE),
  function(x) {
    paste(deparse(formula(x)), collapse = "")
  }
)

# save dredging table results:
m1.cand_tab_clean <- m1.cand_tab_clean %>%
  select(
    formula,
    AICc,
    delta
  ) %>%
  
  mutate(
    across(c(AICc, delta),
           round, 2)
  )


m1.cand_tab_clean

#save
write.csv(m1.cand_tab_clean, 
          "Figures/Calf_Pos_Model_Selection_Fixed_Portion_NND.csv")


#~~~~iii. Residual plots----

#update formula for best model
m1.final <- lme(log(NND) ~ 1 + age_class_coarse + state,
              random = ~ 1 | flight_session_ID,
              data = dat.clust,
              method = "REML")

summary(m1.final)

#open tiff device
png(
  filename = "Figures/Calf_Pos_Model_Diagnostics_Full_NND.png",
  width = 12,
  height = 8,
  units = "in",
  res = 300
)
#
# Create and save multipanel figure for diagnostics:
par(mfrow = c(2, 2))  # 2x2 grid for plots


# Panel 1: Q-Q plot of normalized residuals
qqnorm(resid(m1.final, type = "normalized"), 
       main = "", col = "blue")
qqline(resid(m1.final, type = "normalized"))
title(main = "a)", adj = 0)


# Panel 2: Fitted values vs normalized residuals
plot(fitted(m1.final), resid(m1.final, type = "normalized"),
     xlab = "Fitted Values", ylab = "Normalized Residuals", 
     col = "blue")
title(main = "b)", adj = 0)




# Panel 3: Boxplot of residuals by intensity * treat interaction


boxplot(resid(m1.final, type = "normalized") ~ 
          interaction(dat.clust$age_class_coarse,
                      dat.clust$state,
                      sep = "\n"),
       # xaxt = "n",   # suppress automatic x-axis
        xlab = "Intensity * Treatment", 
        ylab = "Normalized Residuals", 
        cex.axis = 0.5)




title(main = "c)", adj = 0)



# Panel 4: Fitted values vs index to inspect autocorrelation

plot(resid(m1.final, type = "normalized"), 
     ylab = "Normalized Residuals", col = "blue")

title(main = "d)", adj = 0)

dev.off()





#~~~~iv. Effect sizes-----
emm_nnd  <- emmeans(m1.final, ~ age_class_coarse)
emm_nnd  <- emmeans(m1.final, ~ state)

summary(m1.final)


emm_nnd  <- emmeans(m1.final, ~ age_class_coarse | state)

emm_nnd<- as.data.frame(summary(emm_nnd, type = "response"))


emm_nnd$NND <- emm_nnd$response

#~~~~v. Plots-----

p.nnd <- ggplot(
  emm_nnd,
  aes(x = state, y = response)
) +
  
  geom_jitter(
    data = dat.clust,
    inherit.aes = FALSE,
    aes(
      x = state,
      y = NND,
      colour = age_class_coarse
    ),
    position = position_jitterdodge(
      dodge.width = 0.5,
      jitter.width = 0.3
    ),
    alpha = 0.7,
    size = 2
  ) +
  
  geom_point(
    aes(group = age_class_coarse),
    position = position_dodge(width = 0.5),
    size = 2,
    colour = "black"
  ) +
  
  geom_errorbar(
    aes(
      ymin = lower.CL,
      ymax = upper.CL,
      group = age_class_coarse
    ),
    position = position_dodge(width = 0.5),
    width = 0.08,
    colour = "black"
  ) +
  
  scale_color_manual(values = c("gray40", "#69b3a2")) +
  theme_classic() +
  labs(
    x = "State",
    y = "mean NND",
    colour = "Age class",
    title = "a)"
  ) +
  theme(legend.position = "none")
p.nnd


# ~3a. File position -----
#~~~~ i.Select correct random structure (REML)----

# Null model (no random effects)
m2.ls <- gls(filep_f_std ~ 1 + age_class_coarse * state,
             data = dat.clust,
             method = "REML")

summary(m2.ls)

# Model with random intercept for date
m2.ran <- lme(filep_f_std ~ 1 + age_class_coarse * state,
              random = ~ 1 | flight_session_ID,
              data = dat.clust,
              method = "REML")


# Model with both random intercept and variance structure
m2.var.ran <- lme(filep_f_std ~ 1 + age_class_coarse * state,
                  data = dat.clust,
                  random = ~ 1 | flight_session_ID,
                  weights = varIdent(form = ~ 1 | age_class_coarse),
                  method = "REML")


# save
m2.rand.tab <- tibble(
  Model = c(
    "No random structure",
    "Random intercept (flight session)",
    "Random intercept + Variance structure"
  ),
  
  AICc = c(
    AICc(m2.ls),
    AICc(m2.ran),
    AICc(m2.var.ran)
  )
) %>%
  
  mutate(
    deltaAICc = AICc - min(AICc)
  ) %>%
  
  arrange(AICc)

m2.rand.tab

#round to 2 decimal points
m2.rand.tab <- m2.rand.tab %>%
  mutate(
    across(c(AICc, deltaAICc),
           round, 2)
  )

write.csv(m1.rand.tab, "Figures/CalfPos_Model_Selection_Rand_Structure_RFile.csv")



#~~~~ ii.Select Fixed effects (ML)----
m2.final<- update(m2.var.ran, method = "ML")

# Run dredge and keep original object
m2.dredge_results <- dredge(m2.final)

# Convert to dataframe
m2.cand_tab_clean <- as.data.frame(m2.dredge_results)

# Add formulas
m2.cand_tab_clean$formula <- sapply(
  get.models(m2.dredge_results, subset = TRUE),
  function(x) {
    paste(deparse(formula(x)), collapse = "")
  }
)

# save dredging table results:
m2.cand_tab_clean <- m2.cand_tab_clean %>%
  select(
    formula,
    AICc,
    delta
  ) %>%
  
  mutate(
    across(c(AICc, delta),
           round, 2)
  )


m2.cand_tab_clean

#save
write.csv(m2.cand_tab_clean, "Figures/Calf_Pos_Model_Selection_Fixed_Portion_RFile.csv")


#~~~~iii. Residual plots----

#update formula for best model
m2.final <- lme(filep_f_std ~ 1 + age_class_coarse,
                data = dat.clust,
                random = ~ 1 | flight_session_ID,
                weights = varIdent(form = ~ 1 | age_class_coarse),
                method = "ML")


summary(m2.final)

#open tiff device
png(
  filename = "Figures/Calf_Pos_Model_Diagnostics_Full_RFile.png",
  width = 12,
  height = 8,
  units = "in",
  res = 300
)
#
# Create and save multipanel figure for diagnostics:
par(mfrow = c(2, 2))  # 2x2 grid for plots


# Panel 1: Q-Q plot of normalized residuals
qqnorm(resid(m2.final, type = "normalized"), 
       main = "", col = "blue")
qqline(resid(m1.final, type = "normalized"))
title(main = "a)", adj = 0)


# Panel 2: Fitted values vs normalized residuals
plot(fitted(m2.final), resid(m2.final, type = "normalized"),
     xlab = "Fitted Values", ylab = "Normalized Residuals", 
     col = "blue")
title(main = "b)", adj = 0)




# Panel 3: Boxplot of residuals by intensity * treat interaction


boxplot(resid(m2.final, type = "normalized") ~ 
          dat.clust$age_class_coarse,
        # xaxt = "n",   # suppress automatic x-axis
        xlab = "Intensity * Treatment", 
        ylab = "Normalized Residuals", 
        cex.axis = 0.5)




title(main = "c)", adj = 0)



# Panel 4: Fitted values vs index to inspect autocorrelation

plot(resid(m2.final, type = "normalized"), 
     ylab = "Normalized Residuals", col = "blue")

title(main = "d)", adj = 0)

dev.off()





#~~~~iv. Effect sizes-----



emm_rfile  <- emmeans(m2.final, ~ age_class_coarse)

emm_rfile<- as.data.frame(summary(emm_rfile, type = "response"))



#~~~~v. Plots-----


p.f <- ggplot(
  emm_rfile,
  aes(x = age_class_coarse, 
      y = emmean)
  ) +
  
  geom_jitter(
    data = dat.clust,
    inherit.aes = FALSE,
    aes(
      x = age_class_coarse,
      y = filep_f_std,
      colour = age_class_coarse
    ),
    position = position_jitterdodge(
      dodge.width = 0.5,
      jitter.width = 0.3
    ),
    alpha = 0.7,
    size = 2
  ) +
  
  geom_point(
    position = position_dodge(width = 0.5),
    size = 2,
    colour = "black"
  ) +
  
  geom_errorbar(
    aes(
      ymin = lower.CL,
      ymax = upper.CL
    ),
    position = position_dodge(width = 0.5),
    width = 0.02,
    colour = "black"
  ) +
  
  scale_color_manual(values = c("gray40", "#69b3a2")) +
  theme_classic() +
  labs(
    x = "Age class",
    y = "mean RFile",
    colour = "Age class",
    title = "b)"
  ) +
  theme(legend.position = "none")
p.f

# ~3a. Rank position -----
#~~~~ i.Select correct random structure (REML)----

# Null model (no random effects)
m3.ls <- gls(rankP_min_std ~ 1 + age_class_coarse * state,
             data = dat.clust,
             method = "REML")

summary(m3.ls)

# Model with random intercept for date
m3.ran <- lme(rankP_min_std ~ 1 + age_class_coarse * state,
              random = ~ 1 | flight_session_ID,
              data = dat.clust,
              method = "REML")


# Model with both random intercept and variance structure
m3.var.ran <- lme(rankP_min_std ~ 1 + age_class_coarse * state,
                  data = dat.clust,
                  random = ~ 1 | flight_session_ID,
                  weights = varIdent(form = ~ 1 | age_class_coarse),
                  method = "REML")


# save
m3.rand.tab <- tibble(
  Model = c(
    "No random structure",
    "Random intercept (flight session)",
    "Random intercept + Variance structure"
  ),
  
  AICc = c(
    AICc(m3.ls),
    AICc(m3.ran),
    AICc(m3.var.ran)
  )
) %>%
  
  mutate(
    deltaAICc = AICc - min(AICc)
  ) %>%
  
  arrange(AICc)

m3.rand.tab

#round to 2 decimal points
m3.rand.tab <- m3.rand.tab %>%
  mutate(
    across(c(AICc, deltaAICc),
           round, 2)
  )

write.csv(m1.rand.tab, "Figures/CalfPos_Model_Selection_Rand_Structure_RRank.csv")


#~~~~ ii.Select Fixed effects (ML)----
m3.final<- update(m3.var.ran, method = "ML")

# Run dredge and keep original object
m3.dredge_results <- dredge(m3.final)

# Convert to dataframe
m3.cand_tab_clean <- as.data.frame(m3.dredge_results)

# Add formulas
m3.cand_tab_clean$formula <- sapply(
  get.models(m3.dredge_results, subset = TRUE),
  function(x) {
    paste(deparse(formula(x)), collapse = "")
  }
)

# save dredging table results:
m3.cand_tab_clean <- m3.cand_tab_clean %>%
  select(
    formula,
    AICc,
    delta
  ) %>%
  
  mutate(
    across(c(AICc, delta),
           round, 2)
  )


m3.cand_tab_clean

#save
#write.csv(m3.cand_tab_clean, "Figures/Calf_Pos_Model_Selection_Fixed_Portion_RRank.csv")


# keep best: 
# Model with both random intercept and variance structure
m3.final <- lme(rankP_min_std ~ 1 + age_class_coarse,
                  data = dat.clust,
                  random = ~ 1 | flight_session_ID,
                  weights = varIdent(form = ~ 1 | age_class_coarse),
                  method = "ML")


summary(m3.final)

#~~~~iii. Residual plots----



#open tiff device
png(
  filename = "Figures/Calf_Pos_Model_Diagnostics_Full_RRank.png",
  width = 12,
  height = 8,
  units = "in",
  res = 300
)
#
# Create and save multipanel figure for diagnostics:
par(mfrow = c(2, 2))  # 2x2 grid for plots


# Panel 1: Q-Q plot of normalized residuals
qqnorm(resid(m3.final, type = "normalized"), 
       main = "", col = "blue")
qqline(resid(m3.final, type = "normalized"))
title(main = "a)", adj = 0)


# Panel 2: Fitted values vs normalized residuals
plot(fitted(m3.final), resid(m3.final, type = "normalized"),
     xlab = "Fitted Values", ylab = "Normalized Residuals", 
     col = "blue")
title(main = "b)", adj = 0)




# Panel 3: Boxplot of residuals by intensity * treat interaction


boxplot(resid(m3.final, type = "normalized") ~ 
          dat.clust$age_class_coarse,
        
        # xaxt = "n",   # suppress automatic x-axis
        xlab = "Intensity * Treatment", 
        ylab = "Normalized Residuals", 
        cex.axis = 0.5)




title(main = "c)", adj = 0)



# Panel 4: Fitted values vs index to inspect autocorrelation

plot(resid(m3.final, type = "normalized"), 
     ylab = "Normalized Residuals", col = "blue")

title(main = "d)", adj = 0)

dev.off()




#~~~~iv. Effect sizes-----


summary(m3.final)
emm_rrank  <- emmeans(m3.final, ~ age_class_coarse)

emm_rrank<- as.data.frame(summary(emm_rrank, type = "response"))

anova(m3.final)

#~~~~v. Plots-----
p.r <- ggplot(
  emm_rrank,
  aes(x = age_class_coarse, 
      y = emmean)
) +
  
  geom_jitter(
    data = dat.clust,
    inherit.aes = FALSE,
    aes(
      x = age_class_coarse,
      y = rankP_min_std,
      colour = age_class_coarse
    ),
    position = position_jitterdodge(
      dodge.width = 0.5,
      jitter.width = 0.3
    ),
    alpha = 0.7,
    size = 2
  ) +
  
  geom_point(
    position = position_dodge(width = 0.5),
    size = 2,
    colour = "black"
  ) +
  
  geom_errorbar(
    aes(
      ymin = lower.CL,
      ymax = upper.CL
    ),
    position = position_dodge(width = 0.5),
    width = 0.03,
    colour = "black"
  ) +
  
  scale_color_manual(values = c("gray40", "#69b3a2")) +
  theme_classic() +
  labs(
    x = "Age class",
    y = "mean RRank",
    colour = "Age class",
    title = "c)"
  ) +
  theme(legend.position = "none")
p.r


# ~3d. Standardized distance -----
#~~~~ i.Select correct random structure (REML)----

# Null model (no random effects)
m4.ls <- gls(dist_std ~ 1 + age_class_coarse * state,
             data = dat.clust,
             method = "REML")

summary(m4.ls)

# Model with random intercept for date
m4.ran <- lme(dist_std ~ 1 + age_class_coarse * state,
              random = ~ 1 | flight_session_ID,
              data = dat.clust,
              method = "REML")


# Model with both random intercept and variance structure
m4.var.ran <- lme(dist_std ~ 1 + age_class_coarse * state,
                  data = dat.clust,
                  random = ~ 1 | flight_session_ID,
                  weights = varIdent(form = ~ 1 | age_class_coarse),
                  method = "REML")



# save
m4.rand.tab <- tibble(
  Model = c(
    "No random structure",
    "Random intercept (flight session)",
    "Random intercept + Variance structure"
  ),
  
  AICc = c(
    AICc(m4.ls),
    AICc(m4.ran),
    AICc(m4.var.ran)
  )
) %>%
  
  mutate(
    deltaAICc = AICc - min(AICc)
  ) %>%
  
  arrange(AICc)

m4.rand.tab

#round to 2 decimal points
m4.rand.tab <- m4.rand.tab %>%
  mutate(
    across(c(AICc, deltaAICc),
           round, 2)
  )

write.csv(m4.rand.tab, "Figures/CalfPos_Model_Selection_Rand_Structure_SDist.csv")



#~~~~ ii.Select Fixed effects (ML)----
m4.final<- update(m4.var.ran, method = "ML")

# Run dredge and keep original object
m4.dredge_results <- dredge(m4.final)

# Convert to dataframe
m4.cand_tab_clean <- as.data.frame(m4.dredge_results)

# Add formulas
m4.cand_tab_clean$formula <- sapply(
  get.models(m4.dredge_results, subset = TRUE),
  function(x) {
    paste(deparse(formula(x)), collapse = "")
  }
)

# save dredging table results:
m4.cand_tab_clean <- m4.cand_tab_clean %>%
  select(
    formula,
    AICc,
    delta
  ) %>%
  
  mutate(
    across(c(AICc, delta),
           round, 3)
  )


m4.cand_tab_clean

#save
write.csv(m4.cand_tab_clean, "Figures/Calf_Pos_Model_Selection_Fixed_Portion_SDist.csv")


# keep best: 
# Model with both random intercept and variance structure
m4.final <- lme( dist_std~ 1 + age_class_coarse + state ,
                data = dat.clust,
                random = ~ 1 | flight_session_ID,
                weights = varIdent(form = ~ 1 | age_class_coarse),
                method = "ML")


summary(m4.final)

#~~~~iii. Residual plots----



#open tiff device
png(
  filename = "Figures/Calf_Pos_Model_Diagnostics_Full_SDist.png",
  width = 12,
  height = 8,
  units = "in",
  res = 300
)
#
# Create and save multipanel figure for diagnostics:
par(mfrow = c(2, 2))  # 2x2 grid for plots


# Panel 1: Q-Q plot of normalized residuals
qqnorm(resid(m4.final, type = "normalized"), 
       main = "", col = "blue")
qqline(resid(m4.final, type = "normalized"))
title(main = "a)", adj = 0)


# Panel 2: Fitted values vs normalized residuals
plot(fitted(m4.final), resid(m4.final, type = "normalized"),
     xlab = "Fitted Values", ylab = "Normalized Residuals", 
     col = "blue")
title(main = "b)", adj = 0)




# Panel 3: Boxplot of residuals by intensity * treat interaction


boxplot(resid(m4.final, type = "normalized") ~ 
          interaction(dat.clust$age_class_coarse,
                      dat.clust$state),
        
        # xaxt = "n",   # suppress automatic x-axis
        xlab = "Intensity * Treatment", 
        ylab = "Normalized Residuals", 
        cex.axis = 0.5)




title(main = "c)", adj = 0)



# Panel 4: Fitted values vs index to inspect autocorrelation

plot(resid(m4.final, type = "normalized"), 
     ylab = "Normalized Residuals", col = "blue")

title(main = "d)", adj = 0)

dev.off()




#~~~~iv. Effect sizes-----


summary(m4.final)
emm_sdist <- emmeans(m4.final, ~ age_class_coarse)

emm_sdist<- as.data.frame(summary(emm_sdist, type = "response"))

anova(m4.final)

#~~~~v. Plots-----
p.d <- ggplot(
  emm_sdist,
  aes(x = age_class_coarse, 
      y = emmean)
) +
  
  geom_jitter(
    data = dat.clust,
    inherit.aes = FALSE,
    aes(
      x = age_class_coarse,
      y = dist_std,
      
      colour = age_class_coarse
    ),
    position = position_jitterdodge(
      dodge.width = 0.5,
      jitter.width = 0.3
    ),
    alpha = 0.7,
    size = 2
  ) +
  
  geom_point(
    position = position_dodge(width = 0.5),
    size = 2,
    colour = "black"
  ) +
  
  geom_errorbar(
    aes(
      ymin = lower.CL,
      ymax = upper.CL
    ),
    position = position_dodge(width = 0.5),
    width = 0.03,
    colour = "black"
  ) +
  
  scale_color_manual(values = c("gray50", "#69b3a2")) +
  theme_classic() +
  labs(
    x = "Age class",
    y = "mean SDist",
    colour = "Age class",
    title = "d)"
  ) +
  theme(legend.position = "none")
p.d

#significant but tiny?
# 4. SAVE GENERAL OUTPUTS------

#~4a. save plots----
p.out <-  p.nnd + p.f + p.r + p.d

p.out
ggsave("Figures/GLMM_effects_calf_position.png", 
       p.out, width = 8, height = 7)


