#7_Calf position analysis-----
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



#data by flight session
dat<- readRDS(
  "Data/Processed_Data/calf_position_individual_data.RData"
)

# 1. Prepare Data----
# summarize by cluster
dat.block <- dat %>%
  group_by(cluster_ID, age_class_coarse)%>%
  summarise(
    block_id = first(block_id),
    date = first(date), 
    frame_ID = first(frame_ID), 
    phase = first(phase),
    flight_session_ID = first(flight_session_ID),
    intensity_sp = first(intensity_sp),
    NND = mean(NND, na.rm = T),
    filep_f_std = mean(filep_f_std, na.rm = T),
    rankP_min_std = mean(rankP_min_std, na.rm = T),
    dist_std = mean(dist_std, na.rm = T),
    n_total = n(),
  )%>%
  group_by(block_id, age_class_coarse) %>%
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
    n_total = n()
  )





# make a good state factor (i.e., during interaction vs. not)

dat <- dat %>%
  mutate(state = 
           factor(ifelse(is.na(intensity_sp), 
                         "Not",
                         ifelse(intensity_sp == "Neutral", 
                                "Not", 
                                "Aversive"))))




dat.block <- dat.block %>%
  mutate(state = 
           factor(ifelse(is.na(intensity_sp), 
                         "Not",
                         ifelse(intensity_sp == "Neutral", 
                                "Not", 
                                "Aversive"))))

#~2b. general summaries----
dat.block %>%
  filter(age_class_coarse == "Adult")%>%
  group_by(state)%>%
  summarize(n = n())%>%
  ungroup()

dat.block$state <- factor(dat.block$state)



x <- dat %>%
  group_by(cluster_ID)%>%
  summarize(n_calf = sum(age_class_coarse == "Calf", na.rm = T), 
            n_adult = sum(age_class_coarse == "Adult", na.rm = T)) 







# 3. Permutations -----
# within each cluster, re-shuffle age class assignments

n_perm <- 999

out_perms <- vector("list", n_perm)


for(i in seq_len(n_perm)) {
  
  out_perms[[i]] <- dat %>%
    group_by(cluster_ID) %>%
    mutate(age_class_coarse = sample(age_class_coarse)) %>%
    ungroup() %>%
    mutate(permutation = i)
}



output <- bind_rows(out_perms)


# 4. Extract summary statistics ------
# by blocking variable! (phase + flight session)

summ_blocks <- output %>% 
  group_by(permutation, age_class_coarse, block_id) %>%
  summarise(
    state = first(state),
    mean_rank = mean(rankP_min_std, na.rm = TRUE),
    mean_file = mean(filep_f_std, na.rm = TRUE),
    mean_NND = mean(NND),
    n_inds = n(), 
    mean_dist_std = mean(dist_std, na.rm = T),
    .groups = "drop"
  ) %>%
  group_by(permutation, block_id) %>%
  mutate(
    rank_diff_Adult = mean_rank - mean_rank[age_class_coarse == "Adult"],
    file_diff_Adult = mean_file - mean_file[age_class_coarse == "Adult"],
    SDist_diff_Adult = mean_dist_std - mean_dist_std[age_class_coarse == "Adult"],
    NND_diff_Adult = mean_NND - mean_NND[age_class_coarse == "Adult"]
  ) %>%
  ungroup()



# average by permutation

# with state
summ_perms_s <- summ_blocks %>%
  group_by(permutation, 
           age_class_coarse, 
           state)%>%
  summarize(state = first(state), 
            mean_rank = mean(mean_rank), 
            mean_file = mean(mean_file), 
            mean_NND = mean(mean_NND), 
            rank_diff_Adult = mean(rank_diff_Adult),
            file_diff_Adult = mean(file_diff_Adult), 
            NND_diff_Adult = mean(NND_diff_Adult), 
            SDist_diff_Adult= mean(SDist_diff_Adult), 
            .groups = "drop")%>%
  filter(age_class_coarse == "Calf") %>%
  group_by(permutation)%>%
  mutate(int_NND = NND_diff_Adult - NND_diff_Adult[state == "Not"], 
         int_rank = rank_diff_Adult - rank_diff_Adult[state == "Not"],
         int_file = file_diff_Adult - file_diff_Adult[state == "Not"],
         int_edge = SDist_diff_Adult - SDist_diff_Adult[state == "Not"]) %>%
  ungroup()%>%
  filter(state == "Aversive")



# without state
summ_perms <- summ_blocks %>%
  group_by(permutation, 
           age_class_coarse)%>%
  summarize(mean_rank = mean(mean_rank), 
            mean_file = mean(mean_file), 
            mean_NND = mean(mean_NND), 
            mean_dist_std = mean(mean_dist_std),
            rank_diff_Adult = mean(rank_diff_Adult),
            file_diff_Adult = mean(file_diff_Adult), 
            NND_diff_Adult = mean(NND_diff_Adult), 
            SDist_diff_Adult = mean(SDist_diff_Adult), 
            .groups = "drop")



# subset calves
summ_perms_c <- summ_perms %>%
  filter(age_class_coarse == "Calf")


#I. NND -----
# ~~~~~ a. overall ----

nnd.o <- mean(dat.block$NND[dat.block$age_class_coarse == "Calf"] - 
                dat.block$NND[dat.block$age_class_coarse == "Adult"])

nnd.o.sd <- sd(dat.block$NND[dat.block$age_class_coarse == "Calf"] - 
                 dat.block$NND[dat.block$age_class_coarse == "Adult"])


nnd.p <- (sum((summ_perms_c$NND_diff_Adult) <= (nnd.o)) )/ 999


nnd.p


nnd.95<- quantile(summ_perms_c$NND_diff_Adult, probs = c(0.025, 0.975))




p.nnd<- ggplot(data = summ_perms_c, 
       aes(x = NND_diff_Adult))+
  geom_histogram(color = "black", fill = "gray", 
                 binwidth = 0.2)+
  geom_vline(xintercept = nnd.o,  
             col =  "#69b3a2",
             linewidth=1)+
  
  geom_vline(xintercept = nnd.95, 
             lty = "dashed")+
  labs(x = "NND (Calf - Adult)", 
       y = "Count", 
       title = "a)")+
  theme_classic()


# ~~~~~ b. interaction-----
# keep only calves

#difference observed

nnd.diff.aversive <- mean(
  (dat.block$NND[dat.block$age_class_coarse == "Calf" & 
                   dat.block$state == "Aversive"]) - 
  (dat.block$NND[dat.block$age_class_coarse == "Adult" & 
                   dat.block$state == "Aversive"])
  )



nnd.diff.not <- mean( (dat.block$NND[dat.block$age_class_coarse == "Calf" & 
                                       dat.block$state == "Not"]) - 
  (dat.block$NND[dat.block$age_class_coarse == "Adult" & dat.block$state == "Not"])
)


nnd.int.o<- nnd.diff.aversive - nnd.diff.not




# random distribution

nnd.int.95<- quantile(summ_perms_s$int_NND, probs = c(0.025, 0.975))


nnd.int.p <- (sum((summ_perms_s$int_NND) <= (nnd.int.o)) )/ 999


nnd.int.p



p.nnd.int<- ggplot(data = summ_perms_s, 
       aes(x = int_NND))+
  geom_histogram(color = "black", fill = "gray", 
                 bins = 30)+
  geom_vline(xintercept = nnd.int.o,  
             col =  "#69b3a2",
             linewidth=1)+
  
  geom_vline(xintercept = nnd.int.95, 
             lty = "dashed")+
  labs(x = "NND (Calf - Adult) Aversive vs. Not", 
       y = "Count", 
       title = "b)")+
  theme_classic()




#I. Rank -----
# ~~~~~ a. overall ----

rank.o <- mean(dat.block$rankP_min_std[dat.block$age_class_coarse == "Calf"] - 
                dat.block$rankP_min_std[dat.block$age_class_coarse == "Adult"])


rank.o.sd <- sd(dat.block$rankP_min_std[dat.block$age_class_coarse == "Calf"] - 
                 dat.block$rankP_min_std[dat.block$age_class_coarse == "Adult"])

rank.95 <- quantile(summ_perms_c$rank_diff_Adult, probs = c(0.025, 0.975))


rank.p <- (sum((summ_perms_c$rank_diff_Adult) <= (rank.o)) )/ 999



p.rank<-ggplot(data = summ_perms_c, 
       aes(x = rank_diff_Adult))+
  geom_histogram(color = "black", fill = "gray")+
  geom_vline(xintercept = rank.o,  
             col =  "#69b3a2",
             linewidth=1)+
  
  geom_vline(xintercept = rank.95, 
             lty = "dashed")+
  labs(x = "Rank (Calf - Adult)", 
       y = "Count", 
       title = "b)")+
  theme_classic()


# ~~~~~ b. interaction-----

#difference observed

rank.diff.aversive <- mean(
  (dat.block$rankP_min_std[dat.block$age_class_coarse == "Calf" & 
                   dat.block$state == "Aversive"]) - 
    (dat.block$rankP_min_std[dat.block$age_class_coarse == "Adult" & 
                     dat.block$state == "Aversive"])
)



rank.diff.not <- mean( (dat.block$rankP_min_std[dat.block$age_class_coarse == "Calf" & 
                                       dat.block$state == "Not"]) - 
                        (dat.block$rankP_min_std[dat.block$age_class_coarse == "Adult" & dat.block$state == "Not"])
)


rank.int.o<- rank.diff.aversive - rank.diff.not


rank.int.p <- (sum((summ_perms_s$int_rank) <= (rank.int.o)) )/ 999


# random distribution

rank.int.95<- quantile(summ_perms_s$int_rank, probs = c(0.025, 0.975))


p.rank.int <- ggplot(data = summ_perms_s, 
       aes(x = int_rank))+
  geom_histogram(color = "black", fill = "gray", 
                 bins = 30)+
  geom_vline(xintercept = rank.int.o,  
             col =  "#69b3a2",
             linewidth=1)+
  
  geom_vline(xintercept = rank.int.95, 
             lty = "dashed")+
  labs(x = "Rank (Calf - Adult) Aversive vs. Not", 
       y = "Count", 
       title = "b)")+
  theme_classic()



  
#I. File -----
# ~~~~~ a. overall ----

file.o <- mean(dat.block$filep_f_std[dat.block$age_class_coarse == "Calf"] - 
                 dat.block$filep_f_std[dat.block$age_class_coarse == "Adult"])


file.o.sd <- sd(dat.block$filep_f_std[dat.block$age_class_coarse == "Calf"] - 
                 dat.block$filep_f_std[dat.block$age_class_coarse == "Adult"])


file.95 <- quantile(summ_perms_c$file_diff_Adult, probs = c(0.025, 0.975))


file.p <- sum(summ_perms_c$file_diff_Adult >= file.o)/999


p.file <-ggplot(data = summ_perms_c, 
               aes(x = file_diff_Adult))+
  geom_histogram(color = "black", fill = "gray", 
                 bins = 30)+
  geom_vline(xintercept = file.o,  
             col =  "#69b3a2",
             linewidth=1)+
  
  geom_vline(xintercept = file.95, 
             lty = "dashed")+
  labs(x = "File (Calf - Adult)", 
       y = "Count", 
       title = "c)")+
  theme_classic()

p.file

# ~~~~~ b. interaction-----

#difference observed

file.diff.aversive <- mean(
  (dat.block$filep_f_std[dat.block$age_class_coarse == "Calf" & 
                             dat.block$state == "Aversive"]) - 
    (dat.block$filep_f_std[dat.block$age_class_coarse == "Adult" & 
                               dat.block$state == "Aversive"])
)



file.diff.not <- mean( (dat.block$filep_f_std[dat.block$age_class_coarse == "Calf" & 
                                                  dat.block$state == "Not"]) - 
                         (dat.block$filep_f_std[dat.block$age_class_coarse == "Adult" & dat.block$state == "Not"])
)


file.int.o<- file.diff.aversive - file.diff.not


file.int.p <- sum(summ_perms_s$int_file <= file.int.o)/999


# random distribution

file.int.95<- quantile(summ_perms_s$int_file, probs = c(0.025, 0.975))





p.file.int <- ggplot(data = summ_perms_s, 
       aes(x = int_file))+
  geom_histogram(color = "black", fill = "gray", 
                 bins = 30)+
  geom_vline(xintercept = file.int.o,  
             col =  "#69b3a2",
             linewidth=1)+
  
  geom_vline(xintercept = file.int.95, 
             lty = "dashed")+
  labs(x = "File (Calf - Adult) Aversive vs. Not", 
       y = "Count", 
       title = "c)")+
  theme_classic()



#I. Edge -----
# ~~~~~ a. overall ----

edge.o <- mean(dat.block$dist_std[dat.block$age_class_coarse == "Calf"] - 
                 dat.block$dist_std[dat.block$age_class_coarse == "Adult"] )

edge.o.sd <- sd(dat.block$dist_std[dat.block$age_class_coarse == "Calf"] - 
                 dat.block$dist_std[dat.block$age_class_coarse == "Adult"] )


edge.95 <- quantile(summ_perms_c$SDist_diff_Adult, probs = c(0.025, 0.975))


edge.p <- sum(summ_perms_c$SDist_diff_Adult <= edge.o)/999





p.edge <-ggplot(data = summ_perms_c, 
                aes(x = SDist_diff_Adult))+
  geom_histogram(color = "black", fill = "gray", 
                 bins = 30)+
  geom_vline(xintercept = edge.o,  
             col =  "#69b3a2",
             linewidth=1)+
  
  geom_vline(xintercept = edge.95, 
             lty = "dashed")+
  labs(x = "SDist (Calf - Adult)", 
       y = "Count", 
       title = "d)")+
  theme_classic()

p.edge

# ~~~~~ b. interaction-----

#difference observed

edge.diff.aversive <- mean(
  (dat.block$dist_std[dat.block$age_class_coarse == "Calf" & 
                           dat.block$state == "Aversive"]) - 
    (dat.block$dist_std[dat.block$age_class_coarse == "Adult" & 
                             dat.block$state == "Aversive"])
)



edge.diff.not <- mean( (dat.block$dist_std[dat.block$age_class_coarse == "Calf" & 
                                                dat.block$state == "Not"]) - 
                         (dat.block$dist_std[dat.block$age_class_coarse == "Adult" & dat.block$state == "Not"])
)


edge.int.o<- edge.diff.aversive - edge.diff.not

edge.int.p <- sum(summ_perms_s$int_edge <= edge.int.o)/999

# random distribution

edge.int.95<- quantile(summ_perms_s$int_edge, probs = c(0.025, 0.975))


p.edge.int <- ggplot(data = summ_perms_s, 
                 aes(x = int_edge))+
  geom_histogram(color = "black", fill = "gray", 
                 bins = 30)+
  geom_vline(xintercept = edge.int.o,  
             col =  "#69b3a2",
             linewidth=1)+
  
  geom_vline(xintercept = edge.int.95, 
             lty = "dashed")+
  labs(x = "SDist (Calf - Adult) Aversive vs. Not", 
       y = "Count", 
       title = "d)")+
  theme_classic()


 p.edge.int

 
 #Save plot------

 p.out<- (p.nnd + p.rank)/ (p.file + p.edge)

 
 ggsave("Figures/Calf_Adult_Bootstraps.png",
        p.out, width = 8, height = 7)
 
 
 p.out.int<- (p.nnd.int + p.rank.int)/ (p.file.int + p.edge.int)
 
 ggsave("Figures/Calf_Adult_State_Bootstraps.png",
        p.out.int, width = 8, height = 7)
 