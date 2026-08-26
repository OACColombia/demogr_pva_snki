##########################################################
#Spatial PVA for kites template: 
#written by Rob Fletcher and Josh Cullen
##########################################################

#set working directory where data are located
setwd(choose.dir())


library(MASS)
library(expm)   #for matrix powers to calculate growth rate over time
library(tidyverse)
library(wesanderson)
library(future)
library(furrr)
library(zoo) #for rolling sum to determine strings of years


source('helper functions_rf_new.R')

#################
### Load Data ###
#################

# Survival estimates and Transitions among MUs

##survival_mov_PVA code
surv<- read.csv("data/survival_estimates_PVA_age.droughtXpop.inv.csv") %>% 
  mutate(age_cat = case_when(AGE== 'Adult' ~ 'adult',
                             AGE == 'Juvenile' ~ 'juv'))

disp_ad<- read.csv("data/trans_matrix_adults.csv", header = T)
rownames(disp_ad)<- disp_ad$X
disp_ad<- disp_ad[,-1]

new_order <- c("EAST", "EVER", "KRV", "OKEE", "PP", "SJM")

# Reorder both rows and columns
disp_ad <- disp_ad[new_order, new_order]

disp.SE_ad<- read.csv("data/se.trans_matrix_adults.csv", header = T)
rownames(disp.SE_ad)<- disp.SE_ad$X
disp.SE_ad<- disp.SE_ad[,-1]

disp.SE_ad <- disp.SE_ad[new_order, new_order]

disp_juv<- read.csv("data/trans_matrix_juveniles.csv", header = T)
rownames(disp_juv)<- disp_juv$X
disp_juv<- disp_juv[,-1]

disp_juv <- disp_juv[new_order, new_order]

disp.SE_juv<- read.csv("data/se.trans_matrix_juveniles.csv", header = T)
rownames(disp.SE_juv)<- disp.SE_juv$X
disp.SE_juv<- disp.SE_juv[,-1]

disp.SE_juv <- disp.SE_juv[new_order, new_order]

#disp <- read.csv("PVA/trans_noMov.csv")
#rownames(disp)<- disp$X
#disp<- disp[,-1]

#disp.SE <- read.csv("PVA/trans_noMov_se.csv")
#rownames(disp.SE)<- disp.SE$X
#disp.SE<- disp.SE[,-1]


# Reproductive parameters

##Code_nest_success_young_fledged_PVA code
NestSurv <- read.csv("data/PVA_nest_success_estimates.csv", header = T) %>% 
  dplyr::rename(MU = Pop,
                spei_cat = spei,
                invaded = inv.cat,
                NS = ns)

##Code_nest_success_young_fledged_PVA code
Nfledge <- read.csv("data/PVA_young_fledged_estimates.csv", header = T) %>% 
  dplyr::rename( Yng_cnt.SE = se,
                 MU = Pop,
                 spei_cat = spei,
                 invaded = inv.cat)

#PVA_NA_NM
renesting <- read.csv("data/nest_attempts_PVA.csv", header = T) %>%   #pooled across age
  dplyr::rename(Attempt = mean,
                Attempt.SE = se) 


#breedingprob_spei_PVA code
breedprob <- read.csv("data/PVA_breedingprob_est_age_inv.csv", header = T) %>%
  dplyr::rename(Breed = mean,
                spei_cat = spei,
                Breed.SE = se,
                invaded = inv.cat) %>% 
  mutate(age_cat = case_when(age_cat== 'A' ~ 'adult',
                             age_cat == 'J' ~ 'juv'))


# Load minimum number of sighted (banded + unbanded) birds per year and MU for K
K.annual <- read.csv("data/min.known.alive.forPVA.2Apr2025.csv") %>% 
  dplyr::select(-X) 

K.annual<- K.annual %>% 
  filter(Year >= 2019)   #subset to only include 2019-present so representative of current conditions (i.e., exotic snails)

#detectability ratio: adjust based on difference between pop estimate and MNKA
K.annual[is.na(K.annual)] <- 0
K.annual$det <- round(rowSums(K.annual[,2:7], na.rm = T)/K.annual$PopSize, 3)  
K.annual$det <- ifelse(K.annual$det > 1, 1, K.annual$det) #one year slightly > 1
K.annual[,2:7] <- round(K.annual[,2:7] / K.annual$det, 0)  

#proportion for starting distributions
K.annual.prop <- K.annual
K.annual.prop[is.na(K.annual)] <- 0
K.annual.prop <- cbind(K.annual$Year, K.annual.prop[, 2:7] / rowSums(K.annual.prop[, 2:7]))



#####################
### Reformat Data ###
#####################

# Create DF for demographic rates
survNo <- surv %>% 
  filter(drought == "No")

survNormal <- survNo %>% 
  mutate (spei_cat = case_when (drought == "No" ~ "Normal"))

survWet <- survNo %>% 
  mutate (spei_cat = case_when (drought == "No" ~ "Wet"))

survDrought <- surv %>% 
  filter(drought == "Yes") %>% 
  mutate (spei_cat = case_when (drought == "Yes" ~ "Drought"))

surv2 <- rbind(survDrought, survNormal, survWet)

demo.df<- surv2 %>% 
  dplyr::rename(MU = stratum, S = estimate, S.SE = se, S_LCL = lcl, S_UCL = ucl) %>% 
  dplyr::select(-X,
                -AGE) %>% 
  arrange(MU, spei_cat, age_cat) %>% 
  relocate(age_cat, .after = MU) 



# Add estimates from repro params
demo.df2<- demo.df %>% 
  full_join(., breedprob[,c("MU","spei_cat","age_cat","Breed","Breed.SE", "invaded")],
            by = c("MU","spei_cat","age_cat", "invaded")) %>%  #need to decide what to do with age_cat beacuse they are different form survival
  left_join(., NestSurv[,c("MU","spei_cat","NS","NS_SE","ns_lower","ns_upper", "invaded")],
            by = c("MU","spei_cat", "invaded")) %>% 
  left_join(., Nfledge[,c("MU","spei_cat","Yng_cnt","Yng_cnt.SE", "invaded")],
            by = c("MU","spei_cat", "invaded")) %>% 
  left_join(., renesting[,c("MU","spei_cat","Attempt","Attempt.SE")],
            by = c("MU","spei_cat")) %>% 
  filter(invaded == "Yes") %>% 
  dplyr::select(MU, age_cat, spei_cat, S, S.SE, Breed, Breed.SE, 
                NS, NS_SE, Yng_cnt, Yng_cnt.SE, Attempt, Attempt.SE)

#Add NA to Attempts (conditional on breeding)
demo.df2$Attempt<- ifelse(is.na(demo.df2$Breed) | demo.df2$Breed == 0, NA, demo.df2$Attempt)
demo.df2$Attempt.SE<- ifelse(is.na(demo.df2$Breed) | demo.df2$Breed == 0, NA, demo.df2$Attempt.SE)

# Replace NAs with imputed values
#We have no information on PP and droughts or wet (there were no droughts or wet after 2018 when Snail Kites colonized) so we use the same breeding probability than for normal years
demo.df2$Breed<- ifelse(demo.df2$MU == "PP" & is.na(demo.df2$Breed),
                        demo.df2$Breed[demo.df2$MU == "PP" & demo.df2$spei_cat == "Normal"],
                        demo.df2$Breed)
demo.df2$Breed.SE<- ifelse(demo.df2$MU == "PP" & is.na(demo.df2$Breed.SE),
                           demo.df2$Breed.SE[demo.df2$MU == "PP" & demo.df2$spei_cat == "Normal"],
                           demo.df2$Breed.SE)

#Rest of NA correspond to no individuals breeding so they are 0
demo.df2$Breed<- ifelse(is.na(demo.df2$Breed) | demo.df2$Breed == 0,
                        1e-5,
                        demo.df2$Breed)
demo.df2$Breed.SE<- ifelse(is.na(demo.df2$Breed.SE) | demo.df2$Breed.SE == 0,
                           1e-5,
                           demo.df2$Breed.SE)

#NAs for NS correspond to no nests
demo.df2$NS<- ifelse(is.na(demo.df2$NS) | demo.df2$NS == 0, 1e-5, demo.df2$NS)
demo.df2$NS_SE<- ifelse(is.na(demo.df2$NS_SE) | demo.df2$NS_SE < 1e-5, 1e-5, demo.df2$NS_SE)

#We have no information on PP and droughts or wet (there were no droughts or wet after 2018 when Snail Kites colonized) so we use the same nest attempts than for normal years
demo.df2$Attempt<- ifelse(demo.df2$MU == "PP" & is.na(demo.df2$Attempt),
                          demo.df2$Attempt[demo.df2$MU == "PP" & demo.df2$spei_cat == "Normal"],
                          demo.df2$Attempt)
demo.df2$Attempt.SE<- ifelse(demo.df2$MU == "PP" & is.na(demo.df2$Attempt.SE),
                             demo.df2$Attempt.SE[demo.df2$MU == "PP" & demo.df2$spei_cat == "Normal"],
                             demo.df2$Attempt.SE)

demo.df2$Attempt<- ifelse(is.na(demo.df2$Attempt), 1, demo.df2$Attempt)
demo.df2$Attempt.SE<- ifelse(is.na(demo.df2$Attempt.SE), 0, demo.df2$Attempt.SE)

#NAs for Yng_cnt correspond to no nests
demo.df2$Yng_cnt<- ifelse(is.na(demo.df2$Yng_cnt), 0, demo.df2$Yng_cnt)
demo.df2$Yng_cnt.SE<- ifelse(is.na(demo.df2$Yng_cnt.SE), 0, demo.df2$Yng_cnt.SE)



##############################################
#   DEFINE DEMOGRAPHY AND MOVEMENT RATES     #
##############################################

# Set key features of population
Np <- 6 # Number of patches
Ns <- 2 # Number of stages

#------------------------------------------#
# Build the vec-permutation matrix (eq 9)
#------------------------------------------#

P <- matrix(0, Ns*Np, Ns*Np)

for (i in 1:Ns){
  for (j in 1:Np){
    E <- matrix(0,Ns,Np)
    E[i,j] <- 1
    P <- P + kronecker(E,t(E))   
  }
}



##########################################################################
##########################################################################
#Population Projections over observed time intervals
##########################################################################
##########################################################################

#1000 total reps: Nrep * Nclimaterep: 20*50 = 1000

Nrep <- 20  # number of draws from demographic rates
Nstart2023 <- 1410  # adult females in 2023
#Nstart2024 <- round(2007/2,0) #adult females in 2024

#Number of projected years
t <- 31             # year 1 is start year (i.e., 2024)
Nclimaterep <- 50   #random draws of climate categories

take <- c(0,0.1,0.2,0.3,0.4,0.5)
nTake <- length(take)



###############################
# SPEI and Drought Setup (Baseline Only)
###############################

set.seed(123)
spei.proj <- vector("list", Nclimaterep)
for(i in 1:Nclimaterep){
  spei.proj[[i]] <- data.frame(baseline = sample(c("Drought", "Normal", "Wet"), 
                                                 t, replace = TRUE, 
                                                 prob = c(0.09, 0.84, 0.07))) #This is the frequency of extreme values in the past years (i.e., no change in climate scenario for the future)
  
}  


###############################
# PREPARE STORAGE FOR RESULTS
###############################



sim_results <- vector("list", Nclimaterep)
for(j in 1:Nclimaterep){
  sim_results[[j]] <- vector("list", Nrep)
  for(k in 1:Nrep){
    sim_results[[j]][[k]] <- vector("list", nTake)
    names(sim_results[[j]][[k]]) <- paste0("take_", take)
  }
}

# Create parallel storage for carrying capacity values
Klow_store <- vector("list", Nclimaterep)
for(j in 1:Nclimaterep){
  Klow_store[[j]] <- vector("list", Nrep)
  for(k in 1:Nrep){
    Klow_store[[j]][[k]] <- vector("list", nTake)
    names(Klow_store[[j]][[k]]) <- paste0("take_", take)
  }
}


# storage for the MU “take” at time = 1
take_store <- vector("list", Nclimaterep)
for(j in seq_len(Nclimaterep)) {
  take_store[[j]] <- vector("list", Nrep)
  for(k in seq_len(Nrep)) {
    take_store[[j]][[k]] <- vector("list", nTake)
    names(take_store[[j]][[k]]) <- paste0("take_", take)
  }
}

# Define the names for the new multiplier scenarios:
multiplier_names <- c("optimistic", "status_quo", "pessimistic")


###################################
### Run matrix projection model ###
###################################

#loop through SPEI replicates 
for (j in 1:Nclimaterep) {  
  #loop through replicate demographic draws
  for (k in 1:Nrep) {         
    
    #random draw of starting distribution based on 2019-2023 distributions
    MU.prop<- K.annual.prop[sample(nrow(K.annual.prop), size = 1), 2:7]
    
    #In case we want to use the 2023 distributions
    #MU.prop_2023<- K.annual.prop[5, 2:7]
    
    #EAST, EVER, KRV, OKEE, PP, SJM
    Na <- round(Nstart2023 * MU.prop * 0.8) #80% adults
    Nj <- round(Nstart2023 * MU.prop * 0.2) #20% sub-adults
    N <-  t(rbind(Nj, Na))
    
    # Random draw of carrying capacity for rep (2019-2023)
    year.ind<- sample(x = 2019:2023, size = t, replace = TRUE)
    
    #create matrix of K for MU x t; randomly sample from distribution based on year
    K<- sapply(year.ind, function(x) K.annual[K.annual$Year == x, c(2:7)] / 2)  #divide by 2 to account only for females
    K<- matrix(unlist(K), nrow = 6, ncol = t, byrow = FALSE)  #convert from list to matrix
    K[is.na(K)]<- 0  #change NAs to 0
    
    
    # Loop over each "take" value
    for (n in 1:nTake) {
      # Distribute "take" proportionally among MUs:
      floaters_1 <- round(take[n] * Nstart2023 * MU.prop, 0)
      
      state0 <- N #starting population size is the same for all take scenarios
      state0_vec <- as.vector(t(state0))
      
      # Create extra lists in storage for the multiplier scenarios:
      sim_results[[j]][[k]][[n]] <- list()
      #store carrying capacity
      Klow_store[[j]][[k]][[n]] <- list()
      # store the number removed per MU
      take_store[[j]][[k]][[n]] <- floaters_1
      
      ###########################
      # Now loop over multiplier scenarios
      ###########################
      for (p in seq_along(multiplier_names)) {
        mult_name <- multiplier_names[p]
        # Draw a multiplier from the appropriate range:
        if(mult_name == "optimistic"){
          mult_val <- 1.25
        } else if(mult_name == "status_quo"){
          mult_val <- 1
        } else if(mult_name == "pessimistic"){
          mult_val <- 0.75
        }
        
        # Compute Klow for this scenario:
        Klow <- round(mult_val * K, 0)
        # Preallocate a matrix to store the Klow values used at each time step:
        Klow_used <- matrix(NA, nrow = Np, ncol = t)
        
        
        sim_state_take <- matrix(NA, nrow = Np * Ns, ncol = t)
        sim_state_take[,1] <- state0_vec #starting size at t =1
        
        
        #loop through time
        for (i in 2:t) {         
          
          print(c( j, k, take[n], mult_name, i))
          
          
          #create matrices for breeding/survival and dispersal
          B.t <- getB(x = demo.df2, spei = spei.proj[[j]]$baseline[i]) #function for creating B matrix
          B.t.float <- getB_floater(x = demo.df2, spei = spei.proj[[j]]$baseline[i]) #function for creating B matrix
          
          #dispersal matrix: 
          M.t <- getM2(disp_ad = t(disp_ad), disp_ad.SE = t(disp.SE_ad),
                       disp_juv = t(disp_juv), disp_juv.SE = t(disp.SE_juv))
          
          #M.t <- getM(disp = t(disp), disp.SE = t(disp.SE))
          # Store Klow for this time step:
          Klow_used[, i] <- Klow[, i]
          
          #vec permutation: demography then dispersal by pop vector by *patches*
          A.t <- t(P) %*% M.t %*% P %*% B.t   
          A.t.float <- t(P) %*% M.t %*% P %*% B.t.float
          
          #vec permutation: demography then dispersal by pop vector by *stages*
          #A.t <-  M.t %*% P %*% B.t %*% t(P)  
          
          # Get population at previous time step and reshape into MU x stage matrix
          state_prev <- round(sim_state_take[, i - 1], 0)
          state_prev_mat <- matrix(state_prev, nrow = Np, ncol = Ns, byrow = TRUE)
          total_prev <- rowSums(state_prev_mat)
          
          # floaters from K:
          flo_K <- total_prev - Klow[,i]  #if positive, above carrying capacity->floaters
          flo_K[flo_K < 0 ] <- 0
          
          # floaters from "take" only applies at i=2:
          flo_all <- if(i==2) flo_K + floaters_1 else flo_K
          
          breeders <- as.numeric(total_prev - flo_all)
          breeders[breeders < 0 ] <- 0
          breed_frac <- breeders / total_prev
          breed_frac[is.nan(breed_frac)] <- 0
          
          
          #get number of breeders and floaters
          ns.mat.breeder <- state_prev_mat * breed_frac #applies to juv/adults
          ns.mat.breeder[is.nan(ns.mat.breeder)] <- 0
          ns.mat.breeder.vec <- round(as.vector(t(ns.mat.breeder)),0)
          ns.mat.floater <- state_prev_mat - ns.mat.breeder 
          ns.mat.floater.vec <- round(as.vector(t(ns.mat.floater)),0)
          
          
          #project new pop size based on demographics of floaters/breeders
          N.t.i.breed <- A.t %*%  ns.mat.breeder.vec
          N.t.i.float <- A.t.float %*% ns.mat.floater.vec
          
          #sum
          sim_state_take[,i] <- round(N.t.i.breed + N.t.i.float,0)
          
        }#end for i
        
        # Summarize Results for This "Take" Value
        # Aggregate across stages to get total MU population per time step.
        MU_pop <- matrix(NA, nrow = Np, ncol = t)
        for (i in 1:t) {
          state_i <- sim_state_take[, i]
          state_i_mat <- matrix(state_i, nrow = Np, ncol = Ns, byrow = TRUE)
          MU_pop[, i] <- rowSums(state_i_mat)
        }
        
        # Store the MU x time matrix for this replicate and take value.
        sim_results[[j]][[k]][[n]][[mult_name]] <- MU_pop
        
        # Store the Klow matrix for this replicate and take:
        Klow_store[[j]][[k]][[n]][[mult_name]] <- Klow_used
        
      }#end for p
      
    }#end for n
    
  } #end for k
  
} #end for j


MU_names <- c("EAST", "EVER", "KRV", "OKEE", "PP", "SJM")

# Flatten the nested sim_results list into a tibble.
df_results <- imap_dfr(sim_results, function(climate_elem, climate_index) {
  imap_dfr(climate_elem, function(demo_elem, dem_rep) {
    imap_dfr(demo_elem, function(take_elem, take_name) {
      # Extract the numeric take value from its name ("take_0", etc.)
      take_val <- as.numeric(gsub("take_", "", take_name))
      imap_dfr(take_elem, function(mult_elem, mult_name) {
        mu_names <- rownames(mult_elem)
        if (is.null(mu_names)) {
          mu_names <- MU_names
        }
        # mult_elem is a matrix: rows = MU, columns = time steps.
        n_time <- ncol(mult_elem)
        tibble(
          MU           = rep(mu_names, each = n_time),
          year         = rep(seq_len(n_time), times = nrow(mult_elem)),
          pop          = as.vector(t(mult_elem)),  # flatten row-wise
          climate_rep  = climate_index,
          dem_rep      = dem_rep,
          take         = take_val,
          management   = mult_name
        )
      })
    })
  })
})

write.csv(df_results, "df_results.csv")

# Filter by management
df_statusquo <- df_results %>% 
  filter(management =="status_quo")

df_opt <- df_results %>% 
  filter(management =="optimistic")

df_pes <- df_results %>% 
  filter(management =="pessimistic")


#Plot of population trajectories at no take and status quo management
df_take0 <- df_statusquo %>% 
  filter(take == 0) %>%
  mutate(calendar_year = 2022 + year)  # so year=1 -> 2024, year=31 -> 2054


df_plot <- df_take0 %>%
  group_by(MU, calendar_year) %>%
  mutate(pop2 = pop*2) %>%
  dplyr::summarize(
    mean_pop = mean(pop2, na.rm = TRUE),
    lower_95 = quantile(pop2, 0.025, na.rm = TRUE),
    upper_95 = quantile(pop2, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

df_start <- df_take0  %>%  
  mutate(pop2 = pop*2) %>%
  filter(year == 1) %>%            # or calendar_year == 2024
  group_by(MU) %>%
  dplyr::summarize(start_pop = mean(pop2, na.rm = TRUE), .groups = "drop")


ggplot(df_plot, aes(x = calendar_year)) +
  # Ribbon for the 95% interval
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95),
              fill = "blue", alpha = 0.2) +
  # Mean line
  geom_line(aes(y = mean_pop), color = "black", size = 1) +
  
  # Dashed red line for the MU-specific mean start_pop
  geom_hline(
    data = df_start,
    aes(yintercept = start_pop),
    linetype = "dashed",
    color = "red"
  ) +
  
  facet_wrap(~ MU, scales = "free_y") +
  
  # Labels and theme
  labs(
    x = "Calendar Year",
    y = "Population size",
    #title = "Population Trajectories (Take = 0)"
  ) +
  scale_x_continuous(breaks = seq(2024, 2054, 5)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.5),
    strip.text = element_text(face = "bold")
  )


#with three take values (the maximum and with no take):

df_take <- df_statusquo %>% 
  filter(take %in% c(0, 0.2, 0.5)) %>% 
  mutate(calendar_year = 2023 + (year - 1))  


df_plot <- df_take %>%
  group_by(MU, calendar_year, take) %>%
  mutate(pop2 = pop*2) %>%
  dplyr::summarize(
    mean_pop = mean(pop2, na.rm = TRUE),
    lower_95 = quantile(pop2, 0.025, na.rm = TRUE),
    upper_95 = quantile(pop2, 0.975, na.rm = TRUE),
    .groups = "drop"
  )


df_start <- df_take %>%  
  mutate(pop2 = pop*2)%>%
  filter(year == 1) %>%
  group_by(MU, take) %>%
  dplyr::summarize(start_pop = mean(pop2, na.rm = TRUE), .groups = "drop")

df_start %>% 
  pivot_wider(values_from = start_pop,
              names_from = take)

df_plot <- df_plot %>% 
  mutate(take = factor(take, levels = c("0", "0.2", "0.5"),
                             labels = c("0", "20%", "50%")))

ggplot(df_plot, aes(x = calendar_year, group = factor(take))) +
  # Ribbon for the 95% interval, with fill based on take value
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95, fill = factor(take)), alpha = 0.2) +
  # Mean population line, with color based on take value
  geom_line(aes(y = mean_pop, color = factor(take)), size = 1) +
  # Dashed horizontal line for mean starting population
  geom_hline(data = df_start, aes(yintercept = start_pop), color = "grey60",
             linetype = "dashed") +
  facet_wrap(~ MU, scales = "free_y") +
  labs(
    x = "Calendar Year",
    y = "Population Size",
    #title = "Population Trajectories at each MU",
    fill = "Breeders\n removed",
    color = "Breeders\n removed"
  ) +
  scale_x_continuous(breaks = seq(2023, 2053, 5)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.5),
    strip.text = element_text(face = "bold")
  )


#add population growth /lambda by MU
N.mu.pop<- df_results %>% 
  mutate(calendar_year = 2022 + year) %>%
  group_by(climate_rep,dem_rep, take, management, MU) %>% 
  mutate(lambda = pop / lag(pop)) %>% 
  mutate(lambda3yr = rollapply(lambda, 3, mean, fill=NA, align = "center", partial=T))

#make sure there are no NaN / infinite
N.mu.pop$lambda[is.infinite(N.mu.pop$lambda)] = NA
N.mu.pop$lambda[is.nan(N.mu.pop$lambda)] = NA
N.mu.pop$lambda3yr[is.infinite(N.mu.pop$lambda3yr)] = NA
N.mu.pop$lambda3yr[is.nan(N.mu.pop$lambda3yr)] = NA

#get summary values across time
N.mu.pop.sum<- N.mu.pop %>%
  filter(calendar_year > 2023) %>% #remove burnin
  group_by(climate_rep,dem_rep, take, management, MU) %>%
  summarize(Ntot_2053 = last(pop),
            lambda_mean = mean(lambda, na.rm = T),
            lambda_geommean = exp(mean(log(lambda), na.rm = T)),
            lambda3yr_mean = mean(lambda3yr, na.rm = T),
            lambda3yr_geommean = exp(mean(log(lambda3yr), na.rm = T))
  ) 

# Define palette for MUs
pal<- c(wes_palette('Darjeeling1', n = 5),
        wes_palette('Darjeeling2', n = 2)[2])
pal<- pal[c(1,4,3,2,5,6)]

par(mfrow = c(1, 1))


#plots at no take
N.mu.pop.sum0 <- N.mu.pop.sum %>% 
  filter (take == 0)

N.mu.pop.sum0 <- N.mu.pop.sum0 %>% 
  mutate(management = factor(management, levels = c("pessimistic", "status_quo", "optimistic"),
                             labels = c("Pessimistic", "Status quo", "Optimistic")))


#Plot population at time 2050 by MU for each management scenario
MU_size_fig <- ggplot(data = N.mu.pop.sum0, aes(x=management, y = Ntot_2054*2, fill= MU)) +
  #geom_point(position = position_jitter(width = 0.03), size = 0.8, shape = 21, alpha = 0.5) +
  geom_violin(width = 0.8, alpha = 0.5, color = "grey40", adjust = 1.5, trim = T) +
  geom_boxplot(width = 0.1, alpha = 0.7, color = "grey40", outlier.shape = NA) +
  scale_fill_manual("", values = pal, guide=FALSE) +
  labs(x = "Management scenario", y = "Population size at 2054") +
  #scale_y_continuous(breaks = seq(0, 1000, by = 250)) +
  facet_wrap( ~ MU, ncol = 2) +
  #coord_flip()+
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 0, hjust=0.5),
        strip.text = element_text(size = 12, face = 'bold'),
        strip.background = element_blank())


#growth rate

#3 year growth rate
MU_growth_fig <- ggplot(data = N.mu.pop.sum0, aes(x=management, y = lambda3yr_geommean, fill = MU)) +
  geom_hline(yintercept = 1, linetype = "dashed")+
  #geom_point(position = position_jitter(width = 0.05), size = 0.5, shape = 21, alpha = 0.3) +
  geom_violin(width = 0.8, alpha = 0.5, color = "grey40", adjust = 1.5, trim = T) +
  geom_boxplot(width = 0.1, alpha = 0.7, color = "grey40", outlier.shape = NA) +
  #scale_color_manual("", values = pal, guide=FALSE) +
  scale_fill_manual("", values = pal, guide=FALSE) +
  labs(x = "", y = "3-year population growth") +
  ylim(c(0.85, 1.1)) +
  facet_wrap( ~ MU, ncol = 2) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 0, hjust=0.5),
        strip.text = element_text(size = 12, face = 'bold'),
        strip.background = element_blank())






## Plot superpopulation estimates ##
df_super_byrep <- df_results %>% 
  mutate(calendar_year = 2022 + year) %>% 
  mutate(pop2 = pop * 2) %>%  
  group_by(climate_rep, dem_rep, calendar_year, take, management) %>% 
  summarize(total_pop = sum(pop2, na.rm = TRUE), .groups = "drop")

df_super <- df_super_byrep %>% 
  group_by(calendar_year, take, management) %>%
  summarize(
    mean_pop = mean(total_pop, na.rm = TRUE),
    lower_95 = quantile(total_pop, 0.025, na.rm = TRUE),
    upper_95 = quantile(total_pop, 0.975, na.rm = TRUE),
    .groups = "drop"
  )



  
#extract k values  
df_Klow <- imap_dfr(Klow_store, function(climate_elem, climate_index) {
  imap_dfr(climate_elem, function(demo_elem, dem_rep) {
    imap_dfr(demo_elem, function(take_elem, take_name) {
      # Convert take name "take_0" to numeric value 0, etc.
      take_val <- as.numeric(gsub("take_", "", take_name))
      imap_dfr(take_elem, function(mult_elem, mult_name) {
        # Determine the number of time steps in this matrix.
        n_time <- ncol(mult_elem)
        # Retrieve row names; if missing, use the default MU_names.
        mu_names <- rownames(mult_elem)
        if (is.null(mu_names)) {
          mu_names <- MU_names
        }
        tibble(
          MU           = rep(mu_names, each = n_time),
          year         = rep(seq_len(n_time), times = length(mu_names)),
          Klow_val     = as.vector(t(mult_elem)),  # flatten row-wise
          climate_rep  = climate_index,
          dem_rep      = dem_rep,
          take         = take_val,
          management   = mult_name
        )
      })
    })
  })
})

write.csv(df_Klow, "df_klow2.csv")

df_K_sum <- df_Klow %>%
  filter(year >1) %>% 
  group_by(climate_rep, dem_rep, year, take, management) %>%
  dplyr::summarize(
    total_K = sum(Klow_val*2, na.rm = TRUE), .groups = "drop") 


df_Klow_avg <- df_K_sum %>%
  group_by(take, management) %>% 
  dplyr::summarize(
    mean_Klow = mean(total_K, na.rm = TRUE)
  ) 


#plot at no take
df_super_take0 <- df_super %>% 
  filter(take == 0)

df_start <- df_super_take0 %>% 
  filter(calendar_year == 2023)

df_Klow_avg0 <- df_Klow_avg %>% 
  filter(take == 0)



df_super_take0 <- df_super_take0  %>% 
  mutate(management = factor(management, levels = c("pessimistic", "status_quo", "optimistic"),
                             labels = c("Pessimistic", "Status quo", "Optimistic")))


df_start <- df_start %>%
  mutate(management = factor(management, 
                             levels = c("pessimistic", "status_quo", "optimistic"),
                             labels = c("Pessimistic", "Status quo", "Optimistic")))


df_Klow_avg0 <- df_Klow_avg0 %>%
  mutate(management = factor(management, 
                             levels = c("pessimistic", "status_quo", "optimistic"),
                             labels = c("Pessimistic", "Status quo", "Optimistic")))


ggplot(df_super_take0, aes(x = calendar_year)) +
  # Ribbon for the 95% interval
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95),
              fill = "blue", alpha = 0.2) +
  # Mean line
  geom_line(aes(y = mean_pop), color = "black", size = 1) +
  
  # Dashed red line for the MU-specific mean start_pop
  geom_hline(
    data = df_start,
    aes(yintercept = mean_pop),
    linetype = "dashed",
    color = "red"
  ) +
  # Dashed horizontal line for mean carrying capacity
  geom_hline(data = df_Klow_avg0, aes(yintercept =  mean_Klow), color = "black",
             linetype = "dashed") +
  
  facet_wrap(~ management) +
  
  # Labels and theme
  labs(
    x = "Calendar Year",
    y = "Population size",
    #title = "Population Trajectories (Take = 0)"
  ) +
  scale_x_continuous(breaks = seq(2023, 2053, 5)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.5),
    strip.text = element_text(face = "bold")
  )


#with three take values at management "status quo"

df_super_take <- df_super %>% 
  filter(take %in% c(0, 0.2, 0.5),
        management == "status_quo") 


df_start <- df_super_take %>% 
  filter(calendar_year == 2023,
         take == 0) #starting population size is the same independently of take


df_Klow_avg <- df_Klow_avg %>% 
  filter(take %in% c(0, 0.2, 0.5),
         management == "status_quo")

df_super_take <- df_super_take %>% 
  mutate(take = factor(take, levels = c("0", "0.2", "0.5"),
                       labels = c("0", "20%", "50%")))

ggplot(df_super_take, aes(x = calendar_year, group = factor(take))) +
  # Ribbon for the 95% interval, with fill based on take value
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95, fill = factor(take)), alpha = 0.2) +
  # Mean population line, with color based on take value
  geom_line(aes(y = mean_pop, color = factor(take)), size = 1) +
  # Dashed horizontal line for mean starting population
  geom_hline(data = df_start, aes(yintercept =  mean_pop), color = "grey60",
             linetype = "dashed") +
  # Dashed horizontal line for mean carrying capacity
  #geom_hline(data = df_Klow_avg, aes(yintercept =  mean_Klow), color = "black",
  #           linetype = "dashed") +
  labs(
    x = "Calendar Year",
    y = "Population Size",
    #title = "Population Trajectories RW",
    fill = "Breeders\n removed",
    color = "Breeders\n removed"
  ) +
  scale_x_continuous(breaks = seq(2023, 2053, 5)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.5),
    strip.text = element_text(face = "bold")
  )




#get mean size
df_super %>% 
  filter(calendar_year == 2053) 


#add population growth /lambda
N.mu <- df_super_byrep %>% 
  group_by(climate_rep, dem_rep, take, management) %>% 
  mutate(lambda = total_pop / lag(total_pop)) %>% 
  mutate(lambda3yr = rollapply(lambda, 3, mean, fill=NA, align = "center", partial=T))

unique(N.mu$lambda)
N.mu$lambda[is.infinite(N.mu$lambda)] = NA
N.mu$lambda[is.nan(N.mu$lambda)] = NA
N.mu$lambda3yr[is.infinite(N.mu$lambda3yr)] = NA
N.mu$lambda3yr[is.nan(N.mu$lambda3yr)] = NA

N.mu.sum<- N.mu %>%
  filter(calendar_year > 2023) %>% #remove burnin
  group_by(climate_rep, dem_rep, calendar_year, take, management) %>%
  summarize(Ntot_2053 = last(total_pop),
            lambda_mean = mean(lambda, na.rm = T),
            lambda_geommean = exp(mean(log(lambda), na.rm = T)),
            lambda3yr_mean = mean(lambda3yr, na.rm = T),
            lambda3yr_geommean = exp(mean(log(lambda3yr), na.rm = T))
  ) 

#get mean 3-year growth 
N.mu.sum %>%
  group_by(take, management) %>%
  summarize(lambda_mean = mean(lambda3yr_geommean, na.rm = T),
            lambda_median = median(lambda3yr_geommean, na.rm = T) )

N.mu %>%
  filter(calendar_year == 2050) %>% 
  group_by(take, management) %>%
  summarize(Ntot= mean(total_pop))

#plot at no take
N.mu.sum0 <- N.mu.sum %>% 
  filter(take == 0)

N.mu.sum0 <- N.mu.sum0  %>% 
  mutate(management = factor(management, levels = c("pessimistic", "status_quo", "optimistic"),
                             labels = c("Pessimistic", "Status quo", "Optimistic")))

Sup_growth_fig <- ggplot(data = N.mu.sum0, aes(x=management, y = lambda3yr_geommean)) +
  geom_hline(yintercept = 1, linetype = "dashed")+
  #geom_point(position = position_jitter(width = 0.05), size = 0.5, shape = 21, alpha = 0.3) +
  geom_violin(width = 1, alpha = 0.6, adjust = 1.5, trim = T, fill = "steelblue") +
  geom_boxplot(width = 0.1, alpha = 0.5, outlier.shape = NA, fill = "lightblue") +
  #scale_color_manual("", values = pal, guide=FALSE) +
  #scale_fill_manual("", values = pal, guide=FALSE) +
  labs(x = "Management scenario", y = "3-year population growth") +
  ylim(c(0.85, 1.1)) +
 # facet_wrap( ~ MU, ncol = 2) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 0, hjust=0.5),
        strip.text = element_text(size = 12, face = 'bold'),
        strip.background = element_blank())




#Plot super-population at time 2053 at different management scenarios
size_fig <- ggplot(data = N.mu.sum0, aes(x=management, y = Ntot_2054)) +
  geom_point(position = position_jitter(width = 0.03), size = 0.8, shape = 21, alpha = 0.5, fill = "grey40") +
  geom_violin(width = 1, alpha = 0.6, adjust = 1.5, trim = T, fill = "steelblue") +
  geom_boxplot(width = 0.1, alpha = 0.5, outlier.shape = NA, fill = "lightblue") +
  labs(x = "Management scenario", y = "Population size at 2053") +
  ylim(c(0,4000))+
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)
  )

#plot at different take values, management status quo 

N.mu.sum_status <- N.mu.sum %>% 
  filter(management == "status_quo")

N.mu.sum_status<- N.mu.sum_status%>% 
  mutate(take = factor(take, levels = c("0","0.1", "0.2","0.3","0.4", "0.5"),
                       labels = c("0","10", "20", "30","40","50")))


size_fig2 <- ggplot(data = N.mu.sum_status, aes(x=factor(take), y = Ntot_2053, fill = factor(take))) +
  geom_point(position = position_jitter(width = 0.03), size = 0.8, shape = 21, alpha = 0.5, fill = "grey40") +
  geom_violin(width = 1, alpha = 0.6, adjust = 1.5, trim = T) +
  geom_boxplot(width = 0.1, alpha = 0.5, outlier.shape = NA) +
  scale_fill_brewer(name = "Take value", 
                    type = "seq", 
                    palette = "Reds") +
  labs(x = "Percentage of breeders removed", y = "Population size at 2053") +
  ylim(c(1000,3300))+
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "none"
  )



#Box plots to compare between takes at 1,5 and 10 years by MU (management status quo)
desired_years <- c(2, 6, 11) #because year 1 is the start year (i.e. 2024)

df_desired <- df_statusquo  %>% 
  filter(year %in% desired_years) %>% 
  mutate(pop2 = pop*2)

df_desired<- df_desired %>% 
  mutate(take = factor(take, levels = c("0","0.1", "0.2","0.3","0.4", "0.5"),
                       labels = c("0","10%", "20%", "30%","40%","50%")))


take_plot <- ggplot(df_desired, aes(x = factor(year), y = pop2, fill = factor(take))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  scale_fill_brewer(name = "Breeders\n removed", 
                    type = "seq", 
                    palette = "Reds") +
  scale_x_discrete(
    breaks = c("2", "6", "11"),                             
    labels = c("2" = "1" ,"6" = "5", "11" = "10")  
  ) +
  labs(
    x = "Time (Years)",
    y = "Population Size",
    fill = "Take Value",
  ) +
  facet_wrap( ~ MU, scales = "free_y", ncol = 2) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 0, hjust=0.5),
        strip.text = element_text(size = 12, face = 'bold'),
        strip.background = element_blank())

#get the means
take_MU_mean <- df_desired %>% 
  group_by(MU, year, take) %>% 
  summarize(meanPop = mean(pop2)) %>% 
  pivot_wider(values_from = meanPop,
              names_from = year)

write.csv(take_MU_mean, "take mean by MU 2023.csv")

#at super population level
df_super_byrep2 <- df_results %>%  
  mutate(pop2 = pop * 2) %>%  
  group_by(climate_rep, dem_rep, year, take, management) %>% 
  summarize(total_pop = sum(pop2, na.rm = TRUE), .groups = "drop")


df_desired_super <- df_super_byrep2  %>% 
  filter(management == "status_quo",
    year %in% desired_years)

df_desired_super<- df_desired_super %>% 
  mutate(take = factor(take, levels = c("0","0.1", "0.2","0.3","0.4", "0.5"),
                       labels = c("0","10%", "20%", "30%","40%","50%")))


take_plot_super <- ggplot(df_desired_super, aes(x = factor(year), y = total_pop, fill = factor(take))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  scale_fill_brewer(name = "Breeders\n removed", 
                    type = "seq", 
                    palette = "Reds") +
  scale_x_discrete(
    breaks = c("2", "6", "11"),                             
    labels = c("2" = "1" ,"6" = "5", "11" = "10")  
  ) +
  labs(
    x = "Time (Years)",
    y = "Population Size",
    fill = "Take Value",
  ) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 0, hjust=0.5),
        strip.text = element_text(size = 12, face = 'bold'),
        strip.background = element_blank())




#plots of the difference in population size at different take values

df_take_diff <- df_desired %>%
  mutate(pop2 = pop*2) %>% 
  group_by(MU, year, climate_rep, dem_rep) %>%
  mutate(pop_0 = pop2[take == 0],
         diffN_from0 = pop2 - pop_0)



diff_plot <- ggplot(df_take_diff, aes(x = factor(year), y = diffN_from0, fill = factor(take))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  ylim(-120,80)+
  scale_fill_brewer(name = "Breeders\n removed", 
                   type = "seq", 
                   palette = "Reds") +
  scale_x_discrete(
    breaks = c("2", "6", "11"),                             # which factor levels to label
    labels = c("2" = "1" ,"6" = "5", "11" = "10")  
  )+
  facet_wrap( ~ MU, ncol = 2) +
  labs(
    x = "Time (years)",
    y = "Difference in population size",
    fill = "Take value"
  ) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 0, hjust=0.5),
        strip.text = element_text(size = 12, face = 'bold'),
        strip.background = element_blank())


#at super population level
df_take_diff_super <- df_desired_super %>%
  group_by( year, climate_rep, dem_rep) %>%
  mutate(pop_0 = total_pop[take == 0],
         diffN_from0 = total_pop - pop_0)


diff_plot_super <- ggplot(df_take_diff_super, aes(x = factor(year), y = diffN_from0, fill = factor(take))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  ylim(-400,300)+
  scale_fill_brewer(name = "Breeders\n removed", 
                    type = "seq", 
                    palette = "Reds") +
  scale_x_discrete(
    breaks = c("2", "6", "11"),                             # which factor levels to label
    labels = c("2" = "1" ,"6" = "5", "11" = "10")  
  ) +
  labs(
    x = "Time (years)",
    y = "Difference in population size",
    fill = "Take value"
  ) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 0, hjust=0.5),
        strip.text = element_text(size = 12, face = 'bold'),
        strip.background = element_blank())




#get the mean values
df_take_diff_mean <- df_take_diff %>% 
  group_by(MU, take, year) %>% 
  summarize(meanDiff = mean(diffN_from0))

#Mean initial values of take per MU



df_removed <- purrr::imap_dfr(take_store, function(clim, climate_rep) {
  purrr::imap_dfr(clim, function(demo, dem_rep) {
    purrr::imap_dfr(demo, function(take_vec, take_name) {
      take_val <- as.numeric(sub("take_","", take_name))
      # coerce to plain numeric and reassign MU names
      tv       <- as.numeric(unlist(take_vec))
      names(tv)<- MU_names
      
      # now one row per MU
      tibble::enframe(tv, name = "MU", value = "removed") %>%
        dplyr::mutate(
          climate_rep = climate_rep,
          dem_rep     = dem_rep,
          take        = take_val
        ) %>%
        dplyr::select(climate_rep, dem_rep, take, MU, removed)
    })
  })
})


unique(df_removed$MU)

write.csv(df_removed, "df_removed.csv")

df_take_init <- df_removed  %>%
  group_by(MU, take) %>%
  summarize(mean_init = mean(removed),
         lower_95 = quantile(removed, 0.025, na.rm = TRUE),
         upper_95 = quantile(removed, 0.975, na.rm = TRUE))


write.csv(df_take_init, "take_init_MU_2023.csv")


df_take_init_super <- df_removed  %>%
  group_by(take,climate_rep,dem_rep) %>%
  summarize(total_rem = sum(removed)) %>% 
  group_by(take) %>% 
  summarise(mean_init = mean(total_rem),
            lower_95 = quantile(total_rem, 0.025, na.rm = TRUE),
            upper_95 = quantile(total_rem, 0.975, na.rm = TRUE))


df_start_take <- df_super_take %>% 
  filter(calendar_year == 2023)

#plot lambda at different take values at 1, 5, and 10 years

N.mu.pop_desired <- N.mu.pop %>% 
  filter(management == "status_quo",
         year %in% desired_years)


N.mu.pop_desired <- N.mu.pop_desired %>% 
  mutate(take = factor(take, levels = c("0","0.1", "0.2","0.3","0.4", "0.5"),
                       labels = c("0","10%", "20%", "30%","40%","50%")))


take_plot_lambda <- ggplot(N.mu.pop_desired , aes(x = factor(year), y = lambda, fill = factor(take))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  #geom_violin(width = 1, alpha = 0.6, adjust = 1.5, trim = T) +
  #ylim(0.4,3)+
  scale_fill_brewer(name = "Breeders\n removed", 
                    type = "seq", 
                    palette = "Reds") +
  geom_hline(yintercept = 1, linetype = "dashed")+
  scale_x_discrete(
    breaks = c("2", "6", "11"),                             
    labels = c("2" = "1" ,"6" = "5", "11" = "10")  
  ) +
  labs(
    x = "Time (Years)",
    y = "Local Growth Rate",
    fill = "Take Value",
  ) +
  facet_wrap( ~ MU, scales = "free_y", ncol = 2) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 0, hjust=0.5),
        strip.text = element_text(size = 12, face = 'bold'),
        strip.background = element_blank())


#at super population 

N.mu_desired <- N.mu %>% 
  mutate(year = calendar_year - 2022) %>% 
  filter(management == "status_quo",
         year %in% desired_years)

pop_take_mean <- N.mu_desired  %>% 
  group_by(take, year) %>% 
  summarise(mean_pop = mean (total_pop))

N.mu_desired <- N.mu_desired %>% 
  mutate(take = factor(take, levels = c("0","0.1", "0.2","0.3","0.4", "0.5"),
                       labels = c("0","10%", "20%", "30%","40%","50%")))


take_plot_lambda_super <- ggplot(N.mu_desired , aes(x = factor(year), y = lambda, fill = factor(take))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  #geom_violin(width = 1, alpha = 0.6, adjust = 1.5, trim = T) +
  #ylim(0.9,1.05)+
  scale_fill_brewer(name = "Breeders\n removed", 
                    type = "seq", 
                    palette = "Reds") +
  geom_hline(yintercept = 1, linetype = "dashed")+
  scale_x_discrete(
    breaks = c("2", "6", "11"),                             
    labels = c("2" = "1" ,"6" = "5", "11" = "10")  
  ) +
  labs(
    x = "Time (Years)",
    y = "Growth Rate",
    fill = "Take Value",
  ) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 0, hjust=0.5),
        strip.text = element_text(size = 12, face = 'bold'),
        strip.background = element_blank())


lambda_take_mean <- N.mu_desired %>% 
  group_by(take, year) %>% 
  summarise(mean_lambda=mean(lambda),
            median_lambda =median(lambda))
