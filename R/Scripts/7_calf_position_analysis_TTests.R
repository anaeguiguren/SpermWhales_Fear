#7_Calf position analysis-----
rm(list = ls())
source("R/Scripts/7_calf_position_analysis_GEEs_data_prep.R")
set.seed(12345)


library(dplyr)
library(patchwork)
library(lubridate)
library(ggplot2)
library(tidyr)





# 1. Prepare Data----
# make a good state factor (i.e., during interaction vs. not)

dat <- dat %>%
  mutate(state = 
           factor(ifelse(is.na(intensity_sp), 
                         "Not",
                         ifelse(intensity_sp == "Neutral", 
                                "Not", 
                                "Aversive"))), 
         block_run = factor(
           paste(
             age_class_coarse, block_id, sep = "_"
           )), 
         block_group = factor(
           paste(age_class_coarse, group_ID, sep = "_")
         ))




#average by block run 


dat.block <- dat %>%
  group_by(block_run) %>%
  summarize(date = first(date),
            age_class_coarse = first(age_class_coarse),
            block_id = first(block_id),
            state = first(state),
            NND = mean(NND, na.rm = T),
            filep_f_std = mean(filep_f_std, na.rm = T), 
            rankP_min_std = mean(rankP_min_std, na.rm = T), 
            prop_edge = sum(outer_25)/ n(), 
            n = n())



#2. General Summaries-----
dat.block %>%
  filter(age_class_coarse == "Adult")%>%
  group_by(state)%>%
  summarize(n = n())



#I. NND -----
t.test(log(NND) ~ age_class_coarse, data = dat.block)

lm.nnd <- lm(log(NND)~ age_class_coarse, data = dat.block)

newdat <- data.frame(age_class_coarse = levels(factor(dat.block$age_class_coarse)))
predict(lm.nnd, newdata = newdat, interval = "confidence")

plot(lm.nnd)

confint(lm.nnd)

nnd.plot<- ggplot(
  data = dat.block, 
  aes(
    x = age_class_coarse,
    y = log(NND), 
    colour = age_class_coarse, 
    fill = age_class_coarse
     )
  ) +
  
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7)+
  
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.5,
      dodge.width  = 1
    ),
    alpha = 0.9, 
    size = 3
  ) +
  scale_fill_manual(values=c("grey", "#69b3a2")) +
  scale_color_manual(values = c("grey" , "#69b3a2"))+
  labs(x = "", 
       y = "mean NND",
       title = "b)") +
  theme_classic()+
  theme(legend.position = "")


nnd.plot




#II. file  -----
t.test(filep_f_std ~ age_class_coarse, data = dat.block)

lm.file <- lm(filep_f_std~ age_class_coarse, data = dat.block)

summary(lm.file)

predict(lm.file, newdata = newdat, interval = "confidence")

plot(lm.file)

confint(lm.nnd)

file.plot<- ggplot(
  data = dat.block, 
  aes(
    x = age_class_coarse,
    y = filep_f_std, 
    colour = age_class_coarse, 
    fill = age_class_coarse
  )
) +
  
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7)+
  
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.5,
      dodge.width  = 1
    ),
    alpha = 0.9, 
    size = 3
  ) +
  scale_fill_manual(values=c("grey", "#69b3a2")) +
  scale_color_manual(values = c("grey" , "#69b3a2"))+
  labs(x = "", 
       y = "mean File position",
       title = "b)") +
  theme_classic()+
  theme(legend.position = "")


file.plot





#II. rank  -----
t.test(rankP_min_std ~ age_class_coarse, data = dat.block)

lm.rank <- lm(rankP_min_std~ age_class_coarse, data = dat.block)

summary(lm.rank)

predict(lm.rank, newdata = newdat, interval = "confidence")

plot(lm.rank)


rank.plot<- ggplot(
  data = dat.block, 
  aes(
    x = age_class_coarse,
    y = rankP_min_std, 
    colour = age_class_coarse, 
    fill = age_class_coarse
  )
) +
  
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7)+
  
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.5,
      dodge.width  = 1
    ),
    alpha = 0.9, 
    size = 3
  ) +
  scale_fill_manual(values=c("grey", "#69b3a2")) +
  scale_color_manual(values = c("grey" , "#69b3a2"))+
  labs(x = "", 
       y = "mean Rank position",
       title = "b)") +
  theme_classic()+
  theme(legend.position = "")


rank.plot





#II. edge  -----

t.test(prop_edge ~ age_class_coarse, data = dat.block)

lm.edge <- glm(prop_edge~ age_class_coarse, data = dat.block, family = binomial)

summary(lm.edge)

predict(lm.edge, newdata = newdat, interval = "confidence", 
        type = "response")

plot(lm.edge)


edge.plot<- ggplot(
  data = dat.block, 
  aes(
    x = age_class_coarse,
    y = prop_edge, 
    colour = age_class_coarse, 
    fill = age_class_coarse
  )
) +
  
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7)+
  
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.5,
      dodge.width  = 1
    ),
    alpha = 0.9, 
    size = 3
  ) +
  scale_fill_manual(values=c("grey", "#69b3a2")) +
  scale_color_manual(values = c("grey" , "#69b3a2"))+
  labs(x = "", 
       y = "mean Rank position",
       title = "b)") +
  theme_classic()+
  theme(legend.position = "")


edge.plot

#######BOXPLOTS-----


#~~~~rank -----
r.plot<- ggplot(
  data = dat, 
  aes(
    x = age_class_coarse,
    y = rankP_min_std, 
    colour = age_class_coarse, 
    fill = age_class_coarse
  )
) +
  
  geom_boxplot(
               alpha = 0.2)+
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.5,
      dodge.width  = 1
    ),
    alpha = 0.5, 
    size = 3
  ) +
  scale_fill_manual(values=c("grey", "#69b3a2")) +
  scale_color_manual(values = c("grey" , "#69b3a2"))+
  labs(x = "", 
       y = "Rank Position",
       title = "a)") +
  theme_classic()+
  theme(legend.position = "")

r.plot


#~~~~file -----

f.plot<- ggplot(
  data = dat, 
  aes(
    x = age_class_coarse,
    y = filep_f_std, 
    colour = age_class_coarse, 
    fill = age_class_coarse
  )
) +
  
  geom_boxplot(
    alpha = 0.2)+
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.5,
      dodge.width  = 1
    ),
    alpha = 0.5, 
    size = 3
  ) +
  scale_fill_manual(values=c("grey", "#69b3a2")) +
  scale_color_manual(values = c("grey" , "#69b3a2"))+
  labs(x = "", 
       y = "File Position",
       title = "b)") +
  theme_classic()+
  theme(legend.position = "")


f.plot



#~~~~NND -----

nnd.plot<- ggplot(
  data = dat, 
  aes(
    x = age_class_coarse,
    y = log(NND), 
    colour = age_class_coarse, 
    fill = age_class_coarse
  )
) +
  
  geom_boxplot(
    alpha = 0.2)+
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.5,
      dodge.width  = 1
    ),
    alpha = 0.5, 
    size = 3
  ) +
  scale_fill_manual(values=c("grey", "#69b3a2")) +
  scale_color_manual(values = c("grey" , "#69b3a2"))+
  labs(x = "", 
       y = "NND",
       title = "c)") +
  theme_classic()+
  theme(legend.position = "")


nnd.plot


#~~~~ Outer_Edge -----

#~~~~~~~raw data-----

out.plot<- dat %>%
  mutate(outer_25 = factor(outer_25))%>%
  ggplot(aes(x = age_class_coarse, 
             fill = age_class_coarse, 
             pattern = outer_25))+
  geom_bar_pattern(position = "fill",
           colour = "black",
           pattern_fill = "black",
           pattern_colour = "black",
           pattern_density = 0.1,
           pattern_spacing = 0.03, 
           alpha = 0.5)+
  scale_fill_manual(values = c("grey", "#69b3a2")) +
  scale_pattern_manual(values = c("none", "stripe")) +
  labs(x = "", 
       y = "total proportion of individuals",
       title = "d)") +
  theme_classic()+
  theme(legend.position = "")


pos.raw <- (r.plot + f.plot)/ (nnd.plot+ out.plot)

ggsave("Figures/Calf_Adult_Positions_Raw.png",
       pos.raw, width = 8, height = 8)
  

#~~~~~~~modelled outputs-----
library(emmeans)
library(ggeffects)
library(boot) # inv.logit estimate


newdat <- data.frame(age_class_coarse = levels(factor(dat$age_class_coarse)), 
                     outer_25 = 0)

plot(ggemmeans(gee.edge.ind.age, 
               terms = "age_class_coarse"))

summary(gee.edge.ind.age)
confint(gee.edge.ind.age)

cc <- coef(summary(gee.edge.ind.age))
ci_tabl <- with(as.data.frame(cc), 
                cbind(lwr=Estimate-1.96*Std.err,
                      upr=Estimate+1.96*Std.err))

rownames(ci_tabl) <- rownames(cc)



mm.fit1 <- model.matrix(terms(gee.edge.ind.age), newdat)# model matrix
newdat$outer_gee = mm.fit1 %*% coef(gee.edge.ind.age) # point estimates

#fit1[[1]]$geese$vbeta is the robust GEE-based variance-covariance matrix 
tvar1.gee <- diag(mm.fit1 %*% tcrossprod(gee.edge.ind.age$geese$vbeta, mm.fit1))  


newdat <- data.frame(newdat, tlo.gee = newdat$outer_gee-2*sqrt(tvar1.gee), 
                     thi.gee = newdat$outer_gee+2*sqrt(tvar1.gee))

newdat <- newdat %>%
  mutate(
    outer_gee_resp = inv.logit(outer_gee), 
    outer_gee_lo = inv.logit(tlo.gee), 
    outer_gee_hi = inv.logit(thi.gee)
  )


out.mod<-ggplot(newdat, 
       aes(x = age_class_coarse, 
           y = outer_gee_resp))+
  geom_point(aes(colour = age_class_coarse), 
             size = 5)+
  geom_errorbar(aes(ymin = outer_gee_lo, 
                    ymax = outer_gee_hi, 
                    colour = age_class_coarse), 
                 linewidth = 1, 
                width = 0.1)+

  scale_colour_manual(values = c("grey", "#69b3a2"))+
  theme_classic()+
  labs(x = "", 
       y = "Modelled P(Outer edge)")+
  theme(legend.position = "")+
  coord_flip()


ggsave("Figures/Effect_Age_P_outer_edge.png",
       out.mod, width = 4, height = 2)



