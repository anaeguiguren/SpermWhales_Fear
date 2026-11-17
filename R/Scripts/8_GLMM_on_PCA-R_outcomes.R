#' Model Selection Protocol:
#' 1. Select correct random structure (with full model) (REML based)
#' 2. Use dredging approach for selecting variables (ML based)
#' 3. Run visual inspection of model diagnostics.
rm(list = ls())


# 0. LOAD DATA AND LIBRARIES----
library(nlme)
library(MuMIn)
library(ggplot2)
library(emmeans)
library(dplyr)
library(patchwork)

#obs_filtered <- readRDS("Data/Processed_Data/obs_filtered_med_cent.rds")
#modified to include clusters observed in that day but much before/later
obs_filtered <- readRDS("Data/Processed_Data/obs_unfiltered_med_cent.rds")


obs_filtered$treat <- obs_filtered$interaction_window

obs_filtered$date <- factor(obs_filtered$date.x)


# Reorganize factor levels
# For intensity_sp: Neutral, Aversive_Gm, Aversive_Kw
obs_filtered$intensity_sp <- factor(obs_filtered$behav_state_int,
                                    levels = c(
                                      "Neutral", "Aversive_Gm", "Aversive_Oo"
                                      ))

# For treat: 'more than 2 hours apart', 'up to two hours apart', 'during'
obs_filtered$treat <- factor(obs_filtered$treat,
                             levels = c("more than 2 hours", 
                                        "within 2hrs", 
                                        "during"))





#0a. SUMMARY TABLE----
#make a summary of 95 CI for each event state 
#order factor levels
#short lables
obs_filtered$treat_short <- dplyr::recode(
  obs_filtered$treat, 
  "more than 2 hours" = ">2hrs",
  "within 2hrs" = "<2hrs",
  "during" = "during"
)


obs_filtered$event <- factor(obs_filtered$event)

#make short event name
obs_filtered$event_short <- dplyr::recode(
  obs_filtered$event,
  "Event_10" = "E10b",
  "Event_10_dawn_interaction" = "E10a",
  "Event_3" = "E03",
  "Event_2" = "E02",
  "Event_5" = "E05",
  "Event_6" = "E06",
  "Event_9" = "E09",
  
)


#phase temporal lable
obs_filtered$phase_when <- factor(paste0(obs_filtered$phase, 
                                  "_", obs_filtered$treat_short))


# sort
obs_filtered$phase_when <- factor(obs_filtered$phase_when, 
                                  levels = c(
                                    "before_>2hrs", 
                                    "before_<2hrs",
                                    "during_during",
                                    "after_<2hrs",
                                    "after_>2hrs"
                                  ))


# 1. MODEL PC1 -----
# ~1a. Select correct random structure (REML based) for PC1 (Dim.1)----

# Null model (no random effects)
m1.ls <- gls(Dim.1 ~ 1 + treat * intensity_sp,
             data = obs_filtered,
             method = "REML")

# Model with random intercept for date
m1.ran <- lme(Dim.1 ~ 1 + treat * intensity_sp,
              random = ~ 1 | flight_session_ID,
              data = obs_filtered,
              method = "REML")

# Model with variance structure for intensity
m1.var <- gls(Dim.1 ~ 1 + treat * intensity_sp,
              data = obs_filtered,
              weights = varIdent(form = ~ 1 | intensity),
              method = "REML")

# Model with both random intercept and variance structure
m1.var.ran <- lme(Dim.1 ~ 1 + treat * intensity_sp,
                  data = obs_filtered,
                  random = ~ 1 | flight_session_ID,
                  weights = varIdent(form = ~ 1 | intensity),
                  method = "REML")

# Compare models using AICc
AICc_table <- AICc(m1.ls, m1.ran, m1.var, m1.var.ran)
print(AICc_table)  # var.ran is the best based on AICc


# save
rand_tab <- tibble(
  Model = c(
    "No random structure",
    "Random intercept (date)",
    "Variance structure (intensity)",
    "Random intercept + Variance structure"
  ),

  
  AICc = c(
    AICc(m1.ls),
    AICc(m1.ran),
    AICc(m1.var),
    AICc(m1.var.ran)
  )
) %>%
  
  mutate(
    deltaAICc = AICc - min(AICc)
  ) %>%
  
  arrange(AICc)

rand_tab

#round to 2 decimal points
rand_tab <- rand_tab %>%
  mutate(
    across(c(AICc, deltaAICc),
           round, 2)
  )

write.csv(rand_tab, "Figures/Model_Selection_Rand_Structure_PC1-R.csv")



# ~1b. Select fixed variables with dredging approach ----
#ML based
m1.final<- update(m1.var.ran, method = "ML")

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
write.csv(m1.cand_tab_clean, "Figures/Model_Selection_Fixed_Portion_PC1-R.csv")


summary(m1.final)




# ~1c. Residual inspection -----

#adjust lables for graphs
obs_filtered$treat_short <- dplyr::recode(
  obs_filtered$treat, 
  "more than 2 hours" = ">2hrs",
  "within 2hrs" = "<2hrs",
  "during" = "during"
)




#adjust lables for graphs
obs_filtered$intensity_short <- dplyr::recode(
  obs_filtered$intensity_sp, 
  "Neutral" = "Neutral",
  "Aversive_Gm" = paste("Aversive", "PW", sep = "\n"),
  "Aversive_Oo" = paste("Aversive", "KW", sep = "\n")
)
# Panel 3: Boxplot of residuals by intensity * treat interaction

labs <- levels(interaction(obs_filtered$treat_short,
                           obs_filtered$intensity_short,
                           sep = "\n"))




#open tiff device
png(
  filename = "Figures/Model_Diagnostics_Full_PC1-R.png",
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
          interaction(obs_filtered$treat_short,
                      obs_filtered$intensity_short,
                      sep = "\n"),
        xaxt = "n",   # suppress automatic x-axis
        xlab = "Intensity * Treatment", 
        ylab = "Normalized Residuals", 
        cex.axis = 0.5)

# draw ALL tick marks
axis(
  side = 1,
  at = seq_along(labs),
  labels = FALSE,
  tck = -0.02
)


# manually draw ALL labels
text(
  x = seq_along(labs),
  y = par("usr")[3] - 0.12* diff(par("usr")[3:4]),
  labels = labs,
  xpd = NA,
  cex = 0.7,
  adj = 0.5
)
title(main = "c)", adj = 0)



# Panel 4: Fitted values vs index to inspect autocorrelation

plot(resid(m1.final, type = "normalized"), 
     ylab = "Normalized Residuals", col = "blue")

title(main = "d)", adj = 0)

dev.off()

# Reset plotting parameters
par(mfrow = c(1, 1))


  
# 2. MODEL PC2 -----
# ~2a. Select correct random structure (REML based) for PC2 (Dim.2)----

# Null model (no random effects)
m2.ls <- gls(Dim.2 ~ 1 + treat * intensity_sp,
             data = obs_filtered,
             method = "REML")

# Model with random intercept for date
m2.ran <- lme(Dim.2 ~ 1 + treat * intensity_sp,
              random = ~ 1 | flight_session_ID,
              data = obs_filtered,
              method = "REML")

# Model with variance structure for intensity
m2.var <- gls(Dim.2 ~ 1 + treat * intensity_sp,
              data = obs_filtered,
              weights = varIdent(form = ~ 1 | intensity),
              method = "REML")

# Model with both random intercept and variance structure
m2.var.ran <- lme(Dim.2 ~ 1 + treat * intensity_sp,
                  data = obs_filtered,
                  random = ~ 1 | flight_session_ID,
                  weights = varIdent(form = ~ 1 | intensity),
                  method = "REML")

# Compare models using AICc
AICc_table <- AICc(m2.ls, m2.ran, m2.var, m2.var.ran)
print(AICc_table)  # var is the best based on AICc


# save
rand_tab <- tibble(
  Model = c(
    "No random structure",
    "Random intercept (date)",
    "Variance structure (intensity)",
    "Random intercept + Variance structure"
  ),

  
  AICc = c(
    AICc(m2.ls),
    AICc(m2.ran),
    AICc(m2.var),
    AICc(m2.var.ran)
  )
) %>%
  
  mutate(
    deltaAICc = AICc - min(AICc)
  ) %>%
  
  arrange(AICc)

rand_tab

#round to 2 decimal points
rand_tab <- rand_tab %>%
  mutate(
    across(c( AICc, deltaAICc),
           round, 2)
  )

write.csv(rand_tab, "Figures/Model_Selection_Rand_Structure_PC2-R.csv")

# ~2b. Select fixed variables with dredging approach ----
#ML based
m2.final <- update(m2.var.ran, method = "ML")

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
write.csv(m2.cand_tab_clean, "Figures/Model_Selection_Fixed_Portion_PC2-R.csv")




# ~2c. Residual inspection -----

#open tiff device
png(
  filename = "Figures/Model_Diagnostics_Full_PC2-R.tiff",
  width = 12,
  height = 8,
  units = "in",
  res = 300
)


# Create and save multipanel figure for diagnostics:
par(mfrow = c(2, 2))  # 2x2 grid for plots


# Panel 1: Q-Q plot of normalized residuals
qqnorm(resid(m2.final, type = "normalized"), 
       main = "", col = "blue")
qqline(resid(m2.final, type = "normalized"))
title(main = "a)", adj = 0)


# Panel 2: Fitted values vs normalized residuals
plot(fitted(m2.final), resid(m2.final, type = "normalized"),
     xlab = "Fitted Values", ylab = "Normalized Residuals", 
     col = "blue")
title(main = "b)", adj = 0)




# Panel 3: Boxplot of residuals by intensity * treat interaction


boxplot(resid(m2.final, type = "normalized") ~ 
          interaction(obs_filtered$treat_short,
                      obs_filtered$intensity_short,
                      sep = "\n"),
        xaxt = "n",   # suppress automatic x-axis
        xlab = "Intensity * Treatment", 
        ylab = "Normalized Residuals", 
        cex.axis = 0.5)

# draw ALL tick marks
axis(
  side = 1,
  at = seq_along(labs),
  labels = FALSE,
  tck = -0.02
)


# manually draw ALL labels
text(
  x = seq_along(labs),
  y = par("usr")[3] - 0.12* diff(par("usr")[3:4]),
  labels = labs,
  xpd = NA,
  cex = 0.7,
  adj = 0.5
)
title(main = "c)", adj = 0)



# Panel 4: Fitted values vs index to inspect autocorrelation

plot(resid(m2.final, type = "normalized"), 
     ylab = "Normalized Residuals", col = "blue")

title(main = "d)", adj = 0)

dev.off()
# Reset plotting parameters
par(mfrow = c(1, 1))


#3.RESIDUAL ACF PLOTS -----

#open tiff device
png(
  filename = "Figures/Residual_ACF_Full_PC1_PC2-R.png",
  width = 8,
  height = 4,
  units = "in",
  res = 300
)

par(mfrow = c(1, 2))
acf(resid(m1.final, type = "normalized"),
    main = "")
title(main = "a)", adj = 0)

acf(resid(m2.final, type = "normalized"),
    main = "")
title(main = "b)", adj = 0)
dev.off()

#4.EMMEANS PLOTS AND TABLES -----

#colour code
risk_cols <- c("#721F81FF", "#CD4071FF",  "#FD9567FF")

#~4a.PC1----
emm.1 <- emmeans(m1.final, ~ treat * intensity_sp)

# Convert to data frame with confidence intervals
emm.1_df <- as.data.frame(emm.1) %>%
  mutate(
    # Create a combined label for plotting
    combination = paste(intensity_sp, treat, sep = " | ")
  )


# make nice well organized lables

# For treat: 'more than 2 hours apart', 'up to two hours apart', 'during'
obs_filtered$treat <- factor(obs_filtered$treat,
                             levels = c("more than 2 hours", 
                                        "within 2hrs", 
                                       "during"))

# give short names
emm.1_df$combination_neat <- dplyr::recode(
  emm.1_df$combination, 
  "Aversive_Gm | more than 2 hours" = "Aversive\n(PW)\n>2hrs",
  "Aversive_Gm | within 2hrs" = "Aversive\n(PW)\n<2hrs",
  "Aversive_Gm | during" = "Aversive\n(PW)\nduring",
  
  "Aversive_Oo | more than 2 hours" = "Aversive\n(KW)\n>2hrs",
  "Aversive_Oo | within 2hrs" = "Aversive\n(KW)\n<2hrs",
  "Aversive_Oo | during" = "Aversive\n(KW)\nduring",
  
  "Neutral | more than 2 hours" = "Neutral\n>2hrs",
  "Neutral | within 2hrs" = "Neutral\n<2hrs",
  "Neutral | during" = "Neutral\nduring",
  
  

)

# order neatly 
emm.1_df$combination_neat <- factor(
  emm.1_df$combination_neat, 
  levels = c(
    "Neutral\n>2hrs",
    "Neutral\n<2hrs",
    "Neutral\nduring",
    
    "Aversive\n(PW)\n>2hrs",
    "Aversive\n(PW)\n<2hrs",
    "Aversive\n(PW)\nduring",
    
    "Aversive\n(KW)\n>2hrs",
    "Aversive\n(KW)\n<2hrs",
    "Aversive\n(KW)\nduring"
  )
)


  


# create forest plot
plot.emm1 <- ggplot(emm.1_df, 
       aes(x = combination_neat, 
           y = emmean,
           color = intensity_sp)) +
  
  # Add horizontal line at zero (or 1 if using odds ratios)
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  
  # Add error bars for 95% CI
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                 height = 0.2, linewidth = 0.8, 
                width = 0.01) +
  
  # Add points for estimates
  geom_point(size = 3) +
  
  # Color scheme matching your previous plots
  scale_color_manual(
    values = c(
      "Neutral" = risk_cols[3],
      "Aversive_Gm" = risk_cols[2],
      "Aversive_Oo" = risk_cols[1]
    )
  ) +
  
  labs(
    x = "",
    y = "PC1-R Estimated marginal mean (95% CI)",
    title = "a)"

  ) +
  
  theme_classic() 




#~4a.PC2----
emm.2 <- emmeans(m2.final, ~ treat * intensity_sp)

# Convert to data frame with confidence intervals
emm.2_df <- as.data.frame(emm.2) %>%
  mutate(
    # Create a combined label for plotting
    combination = paste(intensity_sp, treat, sep = " | ")
  )


# inherit nice well organized lables from pc1
emm.2_df$combination_neat <- emm.1_df$combination_neat

# create forest plot
plot.emm2 <- ggplot(emm.2_df, 
               aes(x = combination_neat, 
                   y = emmean,
                   color = intensity_sp)) +
  
  # Add horizontal line at zero (or 1 if using odds ratios)
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  
  # Add error bars for 95% CI
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                height = 0.2, linewidth = 0.8, 
                width = 0.01) +
  
  # Add points for estimates
  geom_point(size = 3) +
  
  # Color scheme matching your previous plots
  scale_color_manual(
    values = c(
      "Neutral" = risk_cols[3],
      "Aversive_Gm" = risk_cols[2],
      "Aversive_Oo" = risk_cols[1]
    )
  ) +
  
  labs(
    color = "Intensity",
    x = "Intensity * Phase",
    y = "PC2-R Estimated marginal mean (95% CI)",
    title = "b)"
  ) +
  
  theme_classic() 


plot.emm.both<- plot.emm1 / plot.emm2
plot.emm.both

ggsave("Figures/Estimated_Marginal_Means_PC1_PC2-R.png", 
       plot.emm.both,
       width =8, height = 7)



#4. SAVE ALL FOR PLOTTING-----
save.image(file = "Data/Processed_Data/8_GLMM_on_PCA-R_outcomes_results.RData")

