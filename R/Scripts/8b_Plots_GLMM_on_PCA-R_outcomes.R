# LMM model plots


load("Data/Processed_Data/8_GLMM_on_PCA-R_outcomes_results.RData")


#1. GIVE NEAT NAMES TO OBSERVATIONS-----

obs_filtered <- obs_filtered %>%
  mutate(combination = paste(intensity_sp,treat, sep = " | "))


# give short names
obs_filtered$combination_neat <- dplyr::recode(
  obs_filtered$combination, 
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
obs_filtered$combination_neat <- factor(
  obs_filtered$combination_neat, 
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



panel_labels <- c(
  "Neutral" = "Neutral",
  "Aversive_Gm" = "Aversive (GM)",
  "Aversive_Oo" = "Aversive (KW)"
)



# PC1 plot----
p1 <- ggplot(emm.1_df, 
       aes(x = combination_neat, 
                        y = emmean)) +
  
  #add jittered points of actual data
  
  geom_jitter(
    data = obs_filtered, 
    inherit.aes = F,
    aes(x = combination_neat, 
        y = Dim.1, 
        colour = combination_neat), 
    position = position_jitterdodge(
      dodge.width = 0.5,
      jitter.width = 0.4
    ),
    alpha = 0.8,
    size = 2.5
  )+

  # Add horizontal line at zero 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  
  # Add error bars for 95% CI
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                height = 0.2, 
                linewidth = 0.8, 
                width = 0.03) +
  
  # Add points for estimates
  geom_point(size = 3) +
  
  # Color scheme matching your previous plots
  scale_color_manual(
    values = c(
      "Neutral\n>2hrs"   = "gray40", 
      "Neutral\n<2hrs" = "gray40",
      "Neutral\nduring"  = risk_cols[3],
      
      "Aversive\n(PW)\n>2hrs"   = "gray40", 
      "Aversive\n(PW)\n<2hrs" = "gray40",
      "Aversive\n(PW)\nduring" = risk_cols[2],
      
      "Aversive\n(KW)\n>2hrs"   = "gray40", 
      "Aversive\n(KW)\n<2hrs" = "gray40",
      "Aversive\n(KW)\nduring" = risk_cols[1]
    )
  ) +
  facet_wrap(~intensity_sp, 
             nrow = 1,
             labeller = as_labeller(panel_labels),
             scales = "free_x", 
             space = "free_x")+
  
  labs(
    x = "",
    y = "PC1-R Estimated marginal mean (95% CI)",
    title = "a)"
    
  ) +
  
  theme_classic() +
  theme(legend.position = "null")



# PC2 plot----
p2 <- ggplot(emm.2_df, 
             aes(x = combination_neat, 
                 y = emmean)) +
  
  #add jittered points of actual data
  
  geom_jitter(
    data = obs_filtered, 
    inherit.aes = F,
    aes(x = combination_neat, 
        y = Dim.2, 
        colour = combination_neat), 
    position = position_jitterdodge(
      dodge.width = 0.5,
      jitter.width = 0.4
    ),
    alpha = 0.8,
    size = 2.5
  )+
  
  # Add horizontal line at zero 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  
  # Add error bars for 95% CI
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                height = 0.2, 
                linewidth = 0.8, 
                width = 0.03) +
  
  # Add points for estimates
  geom_point(size = 3) +
  
  # Color scheme matching your previous plots
  scale_color_manual(
    values = c(
      "Neutral\n>2hrs"   = "gray40", 
      "Neutral\n<2hrs" = "gray40",
      "Neutral\nduring"  = risk_cols[3],
      
      "Aversive\n(PW)\n>2hrs"   = "gray40", 
      "Aversive\n(PW)\n<2hrs" = "gray40",
      "Aversive\n(PW)\nduring" = risk_cols[2],
      
      "Aversive\n(KW)\n>2hrs"   = "gray40", 
      "Aversive\n(KW)\n<2hrs" = "gray40",
      "Aversive\n(KW)\nduring" = risk_cols[1]
    )
  ) +
  facet_wrap(~intensity_sp, 
             nrow = 1,
             labeller = as_labeller(panel_labels),
             scales = "free_x", 
             space = "free_x")+
  
  labs(
    x = "",
    y = "PC2-R Estimated marginal mean (95% CI)",
    title = "b)"
    
  ) +
  
  theme_classic() +
  theme(legend.position = "null")

p2

p <- p1 / p2
p



ggsave("Figures/Emmeans_observed_PCs_during_not_by_event.png",
       p, width = 9, height = 6.5, 
       dpi = 300)
