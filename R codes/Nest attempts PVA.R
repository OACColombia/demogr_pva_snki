### Number of nest attempts ###
setwd("C:/Users/bios3/OneDrive/Documents/Snail Kite PVA -2024")
library(tidyverse)
library(wesanderson)
library(dplyr)

source('helper functions.R')  #for custom flatten2() function


#################
### Load Data ###
#################

nests <- read.csv("data/breeding_nest_attempt_forPVA_Apr2025.csv")

spei_inv <- read.csv("data/SPEI_per_yr_long_extreme_inv.csv")

site_cat <- read.csv("Josh codes/data/site_attributes3.csv", header = T)

head(nests)
head(site_cat)


#check for incosistencies
nests %>% 
  filter(YOB > resight_year) %>% 
  tally()

#######################
### Data Formatting ###
#######################

#age of bird
nests2 <- nests[nests$age > 0, ]  #remove birds that attempted to nest in birth year
nests2$age_cat <- ifelse(nests2$age > 1, "Ad", "Juv")



#Change sites so they are consistent with site_cat file (to add detectability)
nests3 <- nests2
sort(unique(nests3$resight_loc))
sort(unique(site_cat$Area))

(x1 <- setdiff(nests3$resight_loc,site_cat$Area))
(x2 <- setdiff(site_cat$Area,nests3$resight_loc))
max_len <- max(length(x1), length(x2))
x1 <- c(x1, rep(NA, max_len - length(x1)))
x2 <- c(x2, rep(NA, max_len - length(x2)))
df <- data.frame(x1, x2)
df

df2 <- data.frame(site_cat$Area,site_cat$region,site_cat$NestDect)
df2

site_cat$Area <- gsub('Everglades National Park','Everglades NP',site_cat$Area)
site_cat$Area<- gsub('Loxahatchee National Wildlife Refuge', 'Loxahatchee NWR', site_cat$Area)
nests3$resight_loc <- gsub('Big Cypress National Preserve','Big Cypress',nests3$resight_loc)
nests3$resight_loc <- gsub('STA1E','STA1',nests3$resight_loc)
site_cat$Area <- gsub("St. Johns Marsh","SJM",site_cat$Area)
site_cat$Area <- gsub("Mary A. Mitigation Bank","Mary A Mitigation Bank",site_cat$Area)
site_cat$Area <- gsub("Lake Side Ranch STA","STA Lakeside Ranch",site_cat$Area)
site_cat$Area <- gsub("Rotenberger Wildlife Management Area","Rotenberger WMA",site_cat$Area)
nests3 <- subset(nests3,resight_loc != "North.of.TOHO")
tmp <- c("Lake Cypress",rep(NA,9),0.62)
site_cat <- rbind(site_cat,tmp)
tmp <- c("Kenansville Lake",rep(NA,9),0.62)
site_cat <- rbind(site_cat,tmp)
tmp <- c("Loxahatchee Slough",rep(NA,9),0.49)
site_cat <- rbind(site_cat,tmp)
tmp <- c("C-44 Impoundment",rep(NA,9),0.83)
site_cat <- rbind(site_cat,tmp)
tmp <- c("A-1 FEB",rep(NA,9),0.83)
site_cat <- rbind(site_cat,tmp)
tmp <- c("Lake Hicpochee Impoundment",rep(NA,9),0.83)
site_cat <- rbind(site_cat,tmp)
tmp <- c("Tiger Lake",rep(NA,9),0.62)
site_cat <- rbind(site_cat,tmp)

#merge with MUs and SPEI
head(nests3)
nests.avail <- merge(nests3[,c("nnest", "resight_year", "resight_loc","age_cat","Region", "unique_id")], 
                     site_cat[, c("Area","NestDect")], 
                     by.x = "resight_loc", by.y = "Area")


nests.avail <- merge(nests.avail, spei_inv, 
                     by.x = c("Region", "resight_year"),  by.y = c("Pop", "year"), all.x = T) %>% 
  dplyr::select(-X)

head(nests.avail)

################################################
#Number of nest attempts / individual / year ###
################################################

#check age differences in general

nests.avail %>% 
  group_by(unique_id, age_cat) %>% 
  summarize(Nnest = sum(nnest))%>%  #get number of nests per individual
  filter(Nnest > 0) %>% 
  group_by(age_cat) %>% 
  summarize(Nnest.mean = mean(Nnest)) #get mean


#summarize number of nests by individual*Area*Year

nests.indiv <- nests.avail %>% 
  group_by(unique_id, resight_loc, resight_year, Region, age_cat, spei, inv.cat) %>% #I added the ivnaded category because we want to estimate the nest attempts under invaded and non-invaded scenarios
  summarize(Nnest = sum(nnest)) %>% 
  filter(Nnest >0)


head(nests.indiv)

#adjust for nest detectability per survey using Area (not MU): 
#lake p: 0.62
#pal p: 0.49
#mixed p: 0.50
#STA p: 0.83

#merge detectability
nests.indiv2 <- merge(nests.indiv, 
                      site_cat[, c("Area", "NestDect")], 
                      by.x = "resight_loc", by.y = "Area")

head(nests.indiv2)

nests.indiv2$NestDect <- as.numeric(nests.indiv2$NestDect)

# two opportunities to dectect, hence the square...
nests.indiv2$Nnestp <- nests.indiv2$Nnest / (1 - (1 - nests.indiv2$NestDect)^2) #2 surveys b/c 31-37 days on average a nest is active

nests.indiv2$MU <- nests.indiv2$Region
nests.indiv2$spei_cat <- nests.indiv2$spei

#summarize number of nests by individual*MU (sample statistics)
nests.summary <- nests.indiv2 %>% 
  group_by(MU, spei_cat) %>% 
  summarize(
    Nnest.mean = mean(Nnestp),
    Nnest.SD = sd(Nnestp),
    Nnest.N = length(Nnestp),
    Nnest.SE = Nnest.SD / sqrt(Nnest.N)
  )


## Perform bootstrap resampling to estimate mean and SD

nests.indiv3 <- nests.indiv2 %>% 
  filter(inv.cat == "Yes")

# Need to convert vars to factors so that unobserved MU*SPEI combinations aren't dropped
nests.indiv2<- nests.indiv2 %>% 
  mutate(across(c('MU', 'spei_cat', 'inv.cat'), factor))



# Split data into nested list by MU → spei → inv.cat
nests_tmp <- nests.indiv2 %>%
  split(.$MU) %>%
  map(~ split(.x, .x$spei_cat)) %>%
  map(~ map(.x, ~ split(.x, .x$inv.cat)))

# Flatten with proper naming

tmp2 <- list()

mus <- unique(nests.indiv2$MU)
speis <- unique(nests.indiv2$spei_cat)
invs <- unique(nests.indiv2$inv.cat)

for (mu in mus) {
  for (sp in speis) {
    for (inv in invs) {
      
      # Assign explicit character strings to use in filter
      cur_mu <- as.character(mu)
      cur_spei <- as.character(sp)
      cur_inv <- as.character(inv)
      
      subset_df <- nests.indiv2 %>%
        filter(MU == cur_mu, spei == cur_spei, inv.cat == cur_inv)
      
      key <- paste(cur_mu, cur_spei, cur_inv, sep = "_")
      
      tmp2[[key]] <- subset_df
    }
  }
}


boot.mu<- vector("list", length(tmp2))
names(boot.mu)<- names(tmp2)
samp<- vector()
var<- "Nnestp"
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

# Calculate the number of nests per MU*SPEI combination
Nnest<- nests.indiv2 %>% 
  group_by(MU, spei, inv.cat) %>% 
  tally() %>% 
  ungroup()





boot.summary <- tibble(
  combo = names(boot.mu),
  mean = colMeans(boot.df, na.rm = TRUE),
  se = apply(boot.df, 2, sd, na.rm = TRUE)
) %>%
  separate(combo, into = c("MU", "spei", "inv.cat"), sep = "_") %>%
  complete(MU, spei, inv.cat) %>%  # ensures all combinations are presen
  left_join (Nnest, by = c("MU", "spei", "inv.cat")) %>% 
  mutate(sd = se*sqrt(n))


write.csv(boot.summary, "nest_attempts_PVA.csv")

#plot (pooled; sample stats)
ggplot(nests.summary, aes(x = spei_cat, y = Nnest.mean)) +
  geom_errorbar(aes(ymin = Nnest.mean - Nnest.SE, ymax = Nnest.mean + Nnest.SE), width = 0) +
  geom_point(size = 2) +
  geom_text(aes(label=paste0("N = ", Nnest.N)), y = 1.1, color = "grey20", size=3) +
  labs(x = "SPEI", y = "Number of Nesting Attempts (SE)") +
  scale_y_continuous(limits = c(1.05, 2)) +
  facet_wrap(~ MU) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(face = "bold", size=14, color = "black"),
        axis.text  = element_text(size = 12, color = "black"))

# Bootstrapped samples
ggplot(boot.summary, aes(x = spei, y = mean, color =inv.cat)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0) +
  geom_point(size = 2) +
  geom_text(aes(label=paste0("N = ", n)), y = 1.1, color = "grey20", size=3) +
  labs(x = "SPEI", y = "Number of Nesting Attempts (SE)") +
  scale_y_continuous(limits = c(1.05, 2)) +
  facet_wrap(~ MU) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(face = "bold", size=14, color = "black"),
        axis.text  = element_text(size = 12, color = "black"))


boot.summary <- boot.summary %>% 
  mutate(lower_ci = mean - (1.96*se),
         upper_ci = mean + (1.96*se))


ggplot(boot.summary, aes(x = spei, y = mean, color = as.factor(inv.cat))) +
  geom_point(position = position_dodge(width = 0.1), size = 3) + 
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                width = 0.2, 
                position = position_dodge(width = 0.1)) +
  facet_grid(~ MU) +
  labs(x = "SPEI Category", y = "Number of nesting attempts (CI)", color = "Invaded") +
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))



#final plot (only for invaded = Yes)
boot.summary_inv <- boot.summary%>% 
  filter(inv.cat == "Yes")

ggplot(boot.summary_inv, aes(x = spei, y = mean)) +
  geom_point(position = position_dodge(width = 0.1), size = 3) + 
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                width = 0.2, 
                position = position_dodge(width = 0.1)) +
  facet_grid(~ MU) +
  labs(x = "SPEI Category", y = "Number of nesting attempts (CI)") +
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))




#summarize number of nests by individual*MU*age
nests.age.summary <- nests.indiv2 %>% 
  group_by(MU, age_cat, spei_cat, inv.cat) %>% 
  summarize(
    Nnest.mean = mean(Nnestp),
    Nnest.SD = sd(Nnestp),
    Nnest.N = length(Nnestp),
    Nnest.SE = Nnest.SD / sqrt(Nnest.N)
  )




## Perform bootstrap resampling to estimate mean and SD

# Need to convert vars to factors so that unobserved MU*SPEI*Age combinations aren't dropped
nests.indiv3<- nests.indiv2 %>% 
  mutate(across(age_cat, factor))

tmp.age2 <- list()

mus <- unique(nests.indiv3$MU)
speis <- unique(nests.indiv3$spei)
ages <- unique(nests.indiv3$age_cat)
invs <- unique(nests.indiv3$inv.cat)

for (mu in mus) {
  for (sp in speis) {
    for (age in ages) {
      for (inv in invs) {
        
        # Assign explicit character strings to use in filter
        cur_mu <- as.character(mu)
        cur_spei <- as.character(sp)
        cur_age <- as.character(age)
        cur_inv <- as.character(inv)
        
        subset_df <- nests.indiv3 %>%
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
var<- "Nnestp"
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


Nests.age <- nests.indiv3 %>%
  group_by(MU, spei, age_cat,inv.cat) %>%
  summarise(NNest = sum(Nnest), .groups = "drop")


boot.summary.age <- boot.summary.age %>%
  left_join(Nests.age, by = c("MU", "spei", "age_cat", "inv.cat")) %>%
  mutate(sd = se * sqrt(NNest))


#plot (sample stats)
ggplot(nests.age.summary, 
       aes(x = spei_cat, group = age_cat, color = inv.cat)) +
  geom_errorbar(aes(ymin = Nnest.mean - Nnest.SD,
                    ymax = Nnest.mean + Nnest.SD), width = 0,
                position = position_dodge(width = 0.3)) +
  geom_point(aes(y = Nnest.mean), size = 2, position = position_dodge(width = 0.3)) +
  geom_text(data = nests.age.summary %>% 
              filter(age_cat == 'Ad'), aes(label=paste0("N = ", Nnest.N)), y = 0.9,
            color=wes_palette("Cavalcanti1", 2)[1], size=3, fontface = "bold") +
  geom_text(data = nests.age.summary %>% 
              filter(age_cat == 'Juv'), aes(label=paste0("N = ", Nnest.N)), y = 0.9,
            color=wes_palette("Cavalcanti1", 2)[2], size=3, fontface = "bold") +
  scale_color_manual("", values = wes_palette("Cavalcanti1", 2)) +
  labs(x = "SPEI", y = "Number of nesting attempts (SD)") +
  scale_y_continuous(limits = c(0.9, 2.6)) +
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
              filter(age_cat == 'Ad'), aes(label=paste0("N = ", NNest)), y = 1,
            color=wes_palette("Cavalcanti1", 2)[1], size=3, fontface = "bold") +
  geom_text(data = boot.summary.age %>% 
              filter(age_cat == 'Juv'), aes(label=paste0("N = ", NNest)), y = 1,
            color=wes_palette("Cavalcanti1", 2)[2], size=3, fontface = "bold") +
  scale_color_manual("", values = wes_palette("Cavalcanti1", 2)) +
  labs(x = "SPEI", y = "Number of nesting attempts (SE)") +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(1, 1.9)) +
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


boot.summary.age$age_cat <- factor(boot.summary.age$age_cat,
                                   levels = c("Juv", "Ad"), 
                                   labels = c( "Subadult","Adult"))

ggplot(boot.summary.age, aes(x = spei, y = mean, color = as.factor(inv.cat))) +
  geom_point(position = position_dodge(width = 0.1), size = 3) + 
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                width = 0.2, 
                position = position_dodge(width = 0.1)) +
  facet_grid(age_cat ~ MU) +
  labs(x = "SPEI Category", y = "Nest attempts", color = "Invaded") +
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))



#final plot (only for invaded = Yes)
boot.summary.age_inv <- boot.summary.age %>% 
  filter(inv.cat == "Yes")

ggplot(boot.summary.age_inv, aes(x = spei, y = mean)) +
  geom_point(position = position_dodge(width = 0.1), size = 3) + 
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                width = 0.2, 
                position = position_dodge(width = 0.1)) +
  facet_grid(age_cat ~ MU) +
  labs(x = "SPEI Category", y = "Number of nesting attempts (CI)") +
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))





######################
### Export results ###
######################


write.csv(boot.summary.age_inv, "nest_attempts_age.PVA.csv")
