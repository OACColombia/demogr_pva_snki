### Number of nest attempts ###
setwd("C:/Users/bios3/OneDrive/Documents/Snail Kite PVA -2024")
library(tidyverse)
library(wesanderson)
library(dplyr)

source('helper functions.R')  #for custom flatten2() function


#################
### Load Data ###
#################

nests <- read.csv("breeding_nest_attempt_forPVA_Apr2025.csv")
head(nests)

spei_inv <- read.csv("SPEI_per_yr_long_extreme_inv.csv")

#spei <- read.csv("Josh codes/data/18-month SPEI data.csv", header = T)

site_cat <- read.csv("Josh codes/data/site_attributes3.csv", header = T)

head(nests)
head(site_cat)

nests2 <- nests

nests %>% 
  filter(YOB > resight_year) %>% 
  tally()

#######################
### Data Formatting ###
#######################

#age of bird
nests2 <- nests2[nests2$age > 0, ]  #remove birds that attempted to nest in birth year
nests2$age_cat <- ifelse(nests2$age > 1, "Ad", "Juv")

nests3 <- nests2

################################################
#Number of nest attempts / individual / year ###
################################################
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

head(nests.avail)
nests.avail <- merge(nests.avail, spei_inv, 
                     by.x = c("Region", "resight_year"),  by.y = c("Pop", "year"), all.x = T)

#check age differences in general

nests.avail %>% 
  group_by(unique_id, age_cat) %>% 
  summarize(Nnest = sum(nnest))%>%  #get number of nests per individual
  filter(Nnest > 0) %>% 
  group_by(age_cat) %>% 
  summarize(Nnest.mean = mean(Nnest)) #get mean


#summarize number of nests by individual*Area*Year

head(nests.avail)


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
#Here is up to where I modified the code
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
nests.indiv3<- nests.indiv3 %>% 
  mutate(across(c('MU', 'spei_cat'), factor))

#Split data
tmp<- nests.indiv3 %>% 
  split(nests.indiv3$MU) %>% 
  map(., ~split(., .$spei_cat, drop = FALSE))
tmp2<- flatten2(tmp, keep_names = TRUE, sep = '_')

#tmp <- nests.indiv2 %>% 
#  split(.$MU) %>% 
#  map(~split(., .$spei_cat, drop = FALSE)) %>% 
#  map(~map(., ~split(., .$age_cat, drop = FALSE)))
#
#tmp2 <- flatten2(tmp, keep_names = TRUE, sep = '_')

boot.mu<- vector("list", length(tmp2))
names(boot.mu)<- names(tmp2)
samp<- vector()
var<- "Nnestp"
boot<- 5000  #number of bootstrap samples


#resample
set.seed(123)
for (i in 1:length(tmp2)) {
  for (j in 1:boot) {
    dat<- tmp2[[i]]
    samp[j]<- mean(base::sample(dat[[var]], size = nrow(dat), replace = T), na.rm = T)
  }
  
  boot.mu[[i]]<- samp
}

boot.df<- bind_cols(boot.mu)


# Calculate the number of nests per MU*SPEI combination
Nnest<- nests.indiv3 %>% 
  group_by(MU, spei_cat) %>% 
  tally() %>% 
  ungroup()

boot.summary<- data.frame(MU = rep(sort(unique(nests.indiv3$MU)), each = 3),
                          spei_cat = rep(sort(unique(nests.indiv3$spei_cat)), 6),
                          mean = colMeans(boot.df, na.rm = TRUE),
                          se = apply(boot.df, 2, sd, na.rm = TRUE)
) %>% 
  left_join(., Nnest, by = c('MU','spei_cat')) %>% 
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
ggplot(boot.summary, aes(x = spei_cat, y = mean)) +
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






#summarize number of nests by individual*MU*age
nests.age.summary <- nests.indiv3 %>% 
  group_by(MU, age_cat, spei_cat) %>% 
  summarize(
    Nnest.mean = mean(Nnestp),
    Nnest.SD = sd(Nnestp),
    Nnest.N = length(Nnestp),
    Nnest.SE = Nnest.SD / sqrt(Nnest.N)
  )




## Perform bootstrap resampling to estimate mean and SD

# Need to convert vars to factors so that unobserved MU*SPEI*Age combinations aren't dropped
nests.indiv3<- nests.indiv3 %>% 
  mutate(across(age_cat, factor))

#Split data
tmp.age<- nests.indiv3 %>% 
  split(nests.indiv2$MU) %>% 
  map(., ~split(., .$spei_cat, drop = FALSE)) %>% 
  map_depth(., 2, ~split(., .$age_cat, drop = FALSE))
tmp.age2<- flatten2(tmp.age, keep_names = TRUE, sep = '_') %>% 
  flatten2(., keep_names = TRUE, sep = '_')

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


# Calculate the number of nests per MU*SPEI combination
Nnest.age<- nests.indiv3 %>% 
  group_by(MU, spei_cat, age_cat) %>% 
  tally() %>% 
  ungroup()

nests.indiv3 %>% 
  filter(MU == "EVER",
         spei_cat == "Wet")

boot.summary.age<- data.frame(MU = rep(sort(unique(nests.indiv3$MU)), each = 3*2),
                              spei_cat = rep(rep(sort(unique(nests.indiv3$spei_cat)), 6), each = 2),
                              age_cat = c('Ad','Juv'),
                              mean = colMeans(boot.df.age, na.rm = TRUE),
                              se = apply(boot.df.age, 2, sd, na.rm = TRUE)
) %>% 
  left_join(., Nnest.age, by = c('MU','spei_cat','age_cat')) %>% 
  mutate(sd = se*sqrt(n))

#plot (by age; samples stats)
ggplot(nests.age.summary, 
       aes(x = spei_cat, y = Nnest.mean, group = age_cat, color = age_cat)) +
  geom_errorbar(aes(ymin = Nnest.mean - Nnest.SE, ymax = Nnest.mean + Nnest.SE), width= 0,
                position = position_dodge(width = 0.3)) +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_text(data = nests.age.summary %>% 
              filter(age_cat == 'Ad'), aes(label=paste0("N = ", Nnest.N)), y = 1.08,
            color=wes_palette("Cavalcanti1", 2)[1], size=3, fontface = "bold") +
  geom_text(data = nests.age.summary %>% 
              filter(age_cat == 'Juv'), aes(label=paste0("N = ", Nnest.N)), y = 1.02,
            color=wes_palette("Cavalcanti1", 2)[2], size=3, fontface = "bold") +
  scale_color_manual("", values = wes_palette("Cavalcanti1", 2)) +
  labs(x = "SPEI", y = "Number of Nesting Attempts (SE)") +
  scale_y_continuous(limits = c(1, 2)) +
  facet_wrap(~ MU) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(face = "bold", size=14, color = "black"),
        axis.text  = element_text(size = 12, color = "black"))


#plot (by age; bootstrap)
ggplot(boot.summary.age, 
       aes(x = spei_cat, y = mean, group = age_cat, color = age_cat)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width= 0,
                position = position_dodge(width = 0.3)) +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_text(data = boot.summary.age %>% 
              filter(age_cat == 'Ad'), aes(label=paste0("N = ", n)), y = 1.08,
            color=wes_palette("Cavalcanti1", 2)[1], size=3, fontface = "bold") +
  geom_text(data = boot.summary.age %>% 
              filter(age_cat == 'Juv'), aes(label=paste0("N = ", n)), y = 1.02,
            color=wes_palette("Cavalcanti1", 2)[2], size=3, fontface = "bold") +
  scale_color_manual("", values = wes_palette("Cavalcanti1", 2)) +
  labs(x = "SPEI", y = "Number of Nesting Attempts (SE)") +
  scale_y_continuous(limits = c(1, 2)) +
  facet_wrap(~ MU) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(face = "bold", size=14, color = "black"),
        axis.text  = element_text(size = 12, color = "black"))




######################
### Export results ###
######################
