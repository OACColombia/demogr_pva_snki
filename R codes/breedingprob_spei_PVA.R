
### Breeding probability for PVA ###

library(tidyverse)
library(MuMIn)
library(wesanderson)

source('helper functions_rf_new.R')


#################
### Load data ###
#################

resights <- read.csv("Breeding status PVA 2025.csv")



spei_inv <- read.csv("SPEI_per_yr_long_extreme_inv.csv")



#We don´t need these steps since Meghan already classified the sites 
#check site names against against each data frame

#site_cat <- read.csv("data/site_attributes2_new.csv",header = T)
#resights[!(resights$resight_loc %in% unique(site_cat$Area)), 'resight_loc'] %>% 
#  unique()
#site_cat[!(site_cat$Area %in% unique(resights$resight_loc)), 'Area'] %>% 
#  unique()


# Make changes so that both sets match
#resights$resight_loc<- gsub("SJM", "St. Johns Marsh", resights$resight_loc)
#resights$resight_loc<- gsub("Mary A Mitigation Bank", "Mary A. Mitigation Bank", resights$resight_loc)
#resights$resight_loc<- gsub("Loxahatchee NWR", "Loxahatchee National Wildlife Refuge", resights$resight_loc)
#resights$resight_loc<- gsub("Rotenberger WMA", "Rotenberger Wildlife Management Area", resights$resight_loc)
#resights$resight_loc<- gsub("Big Cypress National Preserve", "Big Cypress", resights$resight_loc)
#resights$resight_loc<- gsub("Everglades NP", "Everglades National Park", resights$resight_loc)
#resights$resight_loc<- gsub("STA Lakeside Ranch", "Lake Side Ranch STA", resights$resight_loc)
#resights$resight_loc<- gsub("STA1E", "STA1", resights$resight_loc)




#######################
### Data Formatting ###
#######################

#only standardized surveys
surveys <- 1:6
resights <- resights[resights$survey_num %in% surveys, ]

#only 1996-2020
resights <- resights[resights$year > 1995, ]

#age at resight
resights$age <- resights$year - resights$YR_1
sort(unique(resights$age))  # ranges from 1 - 30
resights$age_cat <- ifelse(resights$age == 1, "J", "A")  #assign age classes

#calculate resights by individual
resights$breeding <- ifelse(resights$breeding_status_observed == "B", 1, 0)

head(resights)

resigresightsresight.indiv <- resights %>% 
  group_by(unique.id, year, resight_loc) %>% 
  summarize(breed = max(breeding),
            Nbreed = sum(breeding),
            Nresight = n())


#merge resights with MUs and SPEI

regions <- resights %>% 
  distinct(resight_loc, .keep_all = TRUE)

resight.indiv2 <- left_join(resigresightsresight.indiv , regions[c("resight_loc","Region")], by = "resight_loc")

resight.indiv3 <- left_join(resight.indiv2, spei_inv, join_by (Region == Pop, year == year)) %>% 
  dplyr::select (- X) %>% 
  rename(MU = Region )

#number of indiv region/era
indivN <- resight.indiv3 %>% 
  group_by(MU, spei) %>% 
  summarize(Nindiv = length(unique.id),
            Nindivbreed = sum(breed))

#adjust for delta from Reichert et al. (2012)
indivN_delta <- resight.indiv3 %>% 
  group_by(MU, spei, Nresight, inv.cat) %>% 
  summarize(Nindiv = length(unique.id),
            Nindivbreed = sum(breed))

#detection probability of breeding (delta)
delta <- 0.37

#from Azuma et al. 1990, Spotted Owl book
indivN_delta$NDelta_prop <- (indivN_delta$Nindivbreed / indivN_delta$Nindiv) / 
  (1 - (1 - delta)^indivN_delta$Nresight)
indivN_delta$NDelta_prop <- ifelse(indivN_delta$NDelta_prop > 1, 1, indivN_delta$NDelta_prop)
indivN_delta$NDelta <- round(indivN_delta$NDelta_prop * indivN_delta$Nindiv, 0)

#summarize
indivN_delta_summary <- indivN_delta %>% 
  filter(MU != "N") %>% 
  group_by(MU, spei, inv.cat) %>% 
  summarize(Nindiv = sum(Nindiv),
            Nindivbreed = sum(Nindivbreed),
            NindivbreedDelta = sum(NDelta),
            Ndelta_prop.est = mean(NDelta_prop, na.rm = T),
            Ndelta_prop.SD = sd(NDelta_prop, na.rm = T))  




## Perform bootstrap resampling to estimate mean and SD

# Need to convert vars to factors so that unobserved MU*SPEI combinations aren't dropped
indivN_delta<- indivN_delta %>%
  filter(MU != "N") %>% 
  ungroup() %>% 
  mutate(across(MU:inv.cat, factor))

#Split data
# Split data into nested list by MU → spei → inv.cat
nested_tmp <- indivN_delta %>%
  split(.$MU) %>%
  map(~ split(.x, .x$spei)) %>%
  map(~ map(.x, ~ split(.x, .x$inv.cat)))

# Flatten with proper naming
tmp2 <- list()
for (mu in names(nested_tmp)) {
  for (spei in names(nested_tmp[[mu]])) {
    for (inv in names(nested_tmp[[mu]][[spei]])) {
      key <- paste(mu, spei, inv, sep = "_")
      tmp2[[key]] <- nested_tmp[[mu]][[spei]][[inv]]
    }
  }
}

tmp2 <- list()

mus <- unique(indivN_delta$MU)
speis <- unique(indivN_delta$spei)
invs <- unique(indivN_delta$inv.cat)

for (mu in mus) {
  for (sp in speis) {
    for (inv in invs) {
      
      # Assign explicit character strings to use in filter
      cur_mu <- as.character(mu)
      cur_spei <- as.character(sp)
      cur_inv <- as.character(inv)
      
      subset_df <- indivN_delta %>%
        filter(MU == cur_mu, spei == cur_spei, inv.cat == cur_inv)
      
      key <- paste(cur_mu, cur_spei, cur_inv, sep = "_")
      
      tmp2[[key]] <- subset_df
    }
  }
}


boot.mu<- vector("list", length(tmp2))
names(boot.mu)<- names(tmp2)
samp<- vector()
var<- "NDelta_prop"
boot<- 5000  #number of bootstrap samples


#resample
set.seed(123)
for (i in seq_along(tmp2)) {
  for (j in 1:boot) {
    dat <- tmp2[[i]]
    samp[j] <- mean(base::sample(dat[[var]], size = nrow(dat), replace = TRUE), na.rm = TRUE)
  }
  boot.mu[[i]] <- samp
}

boot.df<- bind_cols(boot.mu)  #convert from list to DF


boot.summary <- tibble(
  combo = names(boot.mu),
  mean = colMeans(boot.df, na.rm = TRUE),
  se = apply(boot.df, 2, sd, na.rm = TRUE)
) %>%
  separate(combo, into = c("MU", "spei", "inv.cat"), sep = "_") %>%
  complete(MU, spei, inv.cat) # ensures all combinations are present
 
indivN <- indivN_delta %>%
  group_by(MU, spei, inv.cat) %>%
  summarise(Nindiv = sum(Nindiv), .groups = "drop")

boot.summary <- boot.summary %>%
  left_join(indivN, by = c("MU", "spei", "inv.cat")) %>%
  mutate(sd = se * sqrt(Nindiv))

#plot (sample stats)
ggplot(indivN_delta_summary, aes(x = spei, colour = inv.cat)) +
  geom_errorbar(aes(ymin = Ndelta_prop.est - Ndelta_prop.SD,
                    ymax = Ndelta_prop.est + Ndelta_prop.SD),
                width= 0, position = position_dodge(0.4)) +
  geom_point(aes(y = Ndelta_prop.est, fill = inv.cat), size = 2, pch = 21, position = position_dodge(0.4)) +
  geom_text(aes(label=paste0("N = ", Nindiv)), y = -0.07, color = "grey20", size=3) +
  labs(x = "SPEI", y = "Breeding Probability (SD)") +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(-0.1, 1)) +
  facet_wrap(~ MU) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 12, face="bold"),
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, face = "bold", color = "black"),
        axis.text  = element_text(size = 12, colour = "black"))


#plot (bootstrap) NEED TO CHECK
ggplot(boot.summary, aes(x = spei, colour = inv.cat)) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se),
                width= 0, position = position_dodge(0.4)) +
  geom_point(aes(y = mean, fill = inv.cat), size = 2, pch = 21, position = position_dodge(0.4)) +
  geom_text(aes(label=paste0("N = ", Nindiv)), y = -0.07, color = "grey20", size=3) +
  labs(x = "SPEI", y = "Breeding Probability (SE)") +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(-0.1, 1)) +
  facet_wrap(~ MU) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 12, face="bold"),
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, face = "bold", color = "black"),
        axis.text  = element_text(size = 12, colour = "black"))



#######################
### By age category ###
#######################

resight.indiv.age <- resights %>% 
  group_by(unique.id, year, age_cat, resight_loc) %>% 
  summarize(breed = max(breeding),
            Nbreed = sum(breeding),
            Nresight = n())

#merge resights with MUs and SPEI
resight.indiv.age2 <- left_join(resight.indiv.age , regions[c("resight_loc","Region")], by = "resight_loc")

resight.indiv.age3 <- left_join(resight.indiv.age2, spei_inv, join_by (Region == Pop, year == year)) %>% 
  dplyr::select (- X) %>% 
  rename(MU = Region )


#number of indiv region/era
indivN.age <- resight.indiv.age3 %>% 
  group_by(MU, age_cat, spei, inv.cat) %>% 
  summarize(Nindiv = length(unique.id),
            Nindivbreed = sum(breed))

#adjust for delta from Reichert et al. (2012)
indivN_delta_age <- resight.indiv.age3 %>% 
  group_by(MU, spei, age_cat,inv.cat, Nresight) %>% 
  summarize(Nindiv = length(unique.id),
            Nindivbreed = sum(breed))

#from Azuma et al. 1990, Spotted Owl book
indivN_delta_age$NDelta_prop <- (indivN_delta_age$Nindivbreed / indivN_delta_age$Nindiv) / 
  (1 - (1 - delta)^indivN_delta_age$Nresight)
indivN_delta_age$NDelta_prop <- ifelse(indivN_delta_age$NDelta_prop > 1, 1,
                                       indivN_delta_age$NDelta_prop)
indivN_delta_age$NDelta <- round(indivN_delta_age$NDelta_prop * indivN_delta_age$Nindiv, 0)


#summarize
indivN_delta_age_summary <- indivN_delta_age %>% 
  filter(MU != "N") %>% 
  group_by(age_cat, MU, spei, inv.cat) %>% 
  summarize(Nindiv = sum(Nindiv),
            Nindivbreed = sum(Nindivbreed),
            NindivbreedDelta = sum(NDelta),
            Ndelta_prop.est = mean(NDelta_prop, na.rm = T),
            Ndelta_prop.SD = sd(NDelta_prop, na.rm = T))



## Perform bootstrap resampling to estimate mean and SD

# Need to convert vars to factors so that unobserved MU*SPEI*Age combinations aren't dropped
indivN_delta_age<- indivN_delta_age %>% 
  filter(MU != "N") %>% 
  ungroup() %>% 
  mutate(across(MU:inv.cat, factor))


# Split data into nested list by MU → spei → inv.cat
tmp.age2 <- list()

mus <- unique(indivN_delta_age$MU)
speis <- unique(indivN_delta_age$spei)
ages <- unique(indivN_delta_age$age_cat)
invs <- unique(indivN_delta_age$inv.cat)

for (mu in mus) {
  for (sp in speis) {
    for (age in ages) {
      for (inv in invs) {
        
        # Assign explicit character strings to use in filter
        cur_mu <- as.character(mu)
        cur_spei <- as.character(sp)
        cur_age <- as.character(age)
        cur_inv <- as.character(inv)
        
        subset_df <- indivN_delta_age %>%
          filter(MU == cur_mu, spei == cur_spei, age_cat == cur_age, inv.cat == cur_inv)
        
        key <- paste(cur_mu, cur_spei, cur_age, cur_inv, sep = "_")
        
        tmp.age2[[key]] <- subset_df
      }
    }
  }
}

boot.mu.age<- vector("list", length(tmp.age2))
names(boot.mu.age)<- names(tmp.age2)
samp<- vector()
var<- "NDelta_prop"
boot<- 5000  #number of bootstrap samples


#resample
set.seed(123)
for (i in 1:length(tmp.age2)) {
  for (j in 1:boot) {
    dat<- tmp.age2[[i]]
    
    if (nrow(dat) == 1) {
      
      samp[j]<- mean(dat[[var]], na.rm = T)
      
    } else {
      
      samp[j]<- mean(sample(dat[[var]], size = nrow(dat), replace = T), na.rm = T)
      
    }
  }
  
  boot.mu.age[[i]]<- samp
}

boot.df.age<- bind_cols(boot.mu.age)

boot.summary.age <- tibble(
  combo = names(boot.mu.age),
  mean = colMeans(boot.df.age, na.rm = TRUE),
  se = apply(boot.df.age, 2, sd, na.rm = TRUE)
) %>%
  separate(combo, into = c("MU", "spei", "age_cat", "inv.cat"), sep = "_")


indivN.age <- indivN_delta_age %>%
  group_by(MU, spei, age_cat,inv.cat) %>%
  summarise(Nindiv = sum(Nindiv), .groups = "drop")


#EAST-Drought-invaded does not exist. Under non-invasion there were no juveniles and no individuals breeding
#PP-Drought-invaded does not exist and pre invasion there were no snail kites.
#PP-Wet (2018) had no resights (but Meghan´s original table had 1 resigth but no breeding)
#SJM-Drought-invaded doesn´t exist. Under non-invasion there were no juveniles and no individuals breeding
#KRV-Wet-invaded had no juveniles. Under non-invasion there were 2 juveniles non breeding
#EAST-Wet-invaded had no juveniles breeding
#SJM-Wet-invaded had no juveniles breeding

boot.summary.age <- boot.summary.age %>%
  left_join(indivN.age, by = c("MU", "spei", "age_cat", "inv.cat")) %>%
  mutate(sd = se * sqrt(Nindiv))



#plot (sample stats)
ggplot(indivN_delta_age_summary, 
       aes(x = spei, group = age_cat, color = inv.cat)) +
  geom_errorbar(aes(ymin = Ndelta_prop.est - Ndelta_prop.SD,
                    ymax = Ndelta_prop.est + Ndelta_prop.SD), width = 0,
                position = position_dodge(width = 0.3)) +
  geom_point(aes(y = Ndelta_prop.est), size = 2, position = position_dodge(width = 0.3)) +
  geom_text(data = indivN_delta_age_summary %>% 
              filter(age_cat == 'A'), aes(label=paste0("N = ", Nindiv)), y = -0.07,
            color=wes_palette("Cavalcanti1", 2)[1], size=3, fontface = "bold") +
  geom_text(data = indivN_delta_age_summary %>% 
              filter(age_cat == 'J'), aes(label=paste0("N = ", Nindiv)), y = -0.13,
            color=wes_palette("Cavalcanti1", 2)[2], size=3, fontface = "bold") +
  scale_color_manual("", values = wes_palette("Cavalcanti1", 2)) +
  labs(x = "SPEI", y = "Breeding Probability (SD)") +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(-0.15, 1)) +
  facet_grid(age_cat ~ MU) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 12, face="bold"),
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, face = "bold", color = "black"),
        axis.text  = element_text(size = 12, colour = "black"))



#plot (bootstrap)
ggplot(boot.summary.age, 
       aes(x = spei, group = age_cat, color = inv.cat)) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se), width = 0,
                position = position_dodge(width = 0.3)) +
  geom_point(aes(y = mean), size = 2, position = position_dodge(width = 0.3)) +
  geom_text(data = boot.summary.age %>% 
              filter(age_cat == 'A'), aes(label=paste0("N = ", Nindiv)), y = -0.07,
            color=wes_palette("Cavalcanti1", 2)[1], size=3, fontface = "bold") +
  geom_text(data = boot.summary.age %>% 
              filter(age_cat == 'J'), aes(label=paste0("N = ", Nindiv)), y = -0.13,
            color=wes_palette("Cavalcanti1", 2)[2], size=3, fontface = "bold") +
  scale_color_manual("", values = wes_palette("Cavalcanti1", 2)) +
  labs(x = "SPEI", y = "Breeding Probability (SE)") +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(-0.15, 1)) +
  facet_grid(age_cat ~ MU) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 12, face="bold"),
        panel.grid = element_blank(),
        axis.title = element_text(size = 14, face = "bold", color = "black"),
        axis.text  = element_text(size = 12, colour = "black"))


boot.summary.age <- boot.summary.age %>% 
  mutate(lower_ci = mean - (1.96*se),
         upper_ci = mean + (1.96*se))
boot.summary.age$lower_ci0 <- ifelse(boot.summary.age$lower_ci <0, 0,boot.summary.age$upper_ci)
boot.summary.age$upper_ci1 <- ifelse(boot.summary.age$upper_ci >1, 1 ,boot.summary.age$upper_ci)

boot.summary.age$age_cat <- factor(boot.summary.age$age_cat,
                                   levels = c("J", "A"), 
                                   labels = c( "Subadult","Adult"))

ggplot(boot.summary.age, aes(x = spei, y = mean, color = as.factor(inv.cat))) +
  geom_point(position = position_dodge(width = 0.1), size = 3) + 
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                width = 0.2, 
                position = position_dodge(width = 0.1)) +
  facet_grid(age_cat ~ MU) +
  labs(x = "SPEI Category", y = "Breeding Probability", color = "Invaded") +
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))


#final plot (only for invaded = Yes)
boot.summary.age_inv <- boot.summary.age %>% 
  filter(inv.cat == "Yes")

ggplot(boot.summary.age_inv, aes(x = spei, y = mean, color= age_cat)) +
  geom_point(position = position_dodge(width = 0.1), size = 3) + 
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                width = 0.2, 
                position = position_dodge(width = 0.1)) +
  facet_grid( ~ MU) +
  labs(x = "SPEI Category", y = "Breeding Probability (CI)", color = "Age") +
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))


######################
### Export results ###
######################

# write.csv(boot.summary, "PVA_breedingprob_est.csv", row.names = F)
write.csv(boot.summary.age, "PVA_breedingprob_est_age_inv.csv", row.names = F)
