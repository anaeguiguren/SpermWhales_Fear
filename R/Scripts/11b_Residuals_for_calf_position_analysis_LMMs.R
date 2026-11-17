# List your models
models <- list(m1.final, m2.final, m3.final, m4.final)

# Optional panel labels
panel_labels <- letters[1:8]






#open png
png(
  filename = "Figures/Calf_Pos_Model_Diagnostics_All.png",
  width = 7,
  height = 8,
  units = "in",
  res = 300
)
#

# Create 4 x 2 multipanel figure
par(mfrow = c(4, 2),
    mar = c(2, 4, 2, 1))

for(i in seq_along(models)) {
  
  # Extract normalized residuals
  res <- resid(models[[i]], type = "normalized")
  
  # Show x-axis labels only for bottom row (g and h)
  xlab_text <- if(i == 4) "Index" else ""
  
  # ---- Column 1: Q-Q plot ----
  qqnorm(res,
         main = "",
         col = "blue",
         xlab = xlab_text)
  
  qqline(res)
  
  title(main = paste0(panel_labels[(i - 1) * 2 + 1], ")"),
        adj = 0)
  
  
  # ---- Column 2: Residuals vs Index ----
  plot(res,
       ylab = "Normalized Residuals",
       xlab = xlab_text,
       col = "blue",
       main = "")
  
  title(main = paste0(panel_labels[(i - 1) * 2 + 2], ")"),
        adj = 0)
}

dev.off()
