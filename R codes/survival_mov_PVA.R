### Run multistate CJS model using MARK via RMark to calculate survival and transitions among populations 

#survival: can vary by age, population, spei, and invasion of P. maculata
#movement: can vary by age and distance between populations
#detection: varies by wetland type

library(RMark)  
library(dplyr)
library(tidyverse)
library(ggplot2)
library(car)

#see README for column descriptions

# Load capture histories
ms <- read.csv("survival_movement_data_Feb2024.csv")
#to make single vector detection history (for MARK)
ms <- ms[,-1]
ms$ch <- do.call(paste, c(ms[,1:length(ms)], sep=""))
ms$ch <- as.character(ms$ch)
t1<- 1996 # starting year for encounter histories

head(ms)

ms.pr<- process.data(ms, begin.time = t1, model = "Multistrata")
head(ms.pr$data)

# Create design data
ms.ddl<- make.design.data(ms.pr)

#If we decide not to include AGE into the movement model, we can set pim.type="time" to reduce the number of parameters and make the model run much faster (like Laake suggested to speed up and Josh did in the past) 
#ms.ddl<- make.design.data(ms.pr, parameters=list(Psi=list(pim.type="time")))


# Add age classes (juvenile, adult) for survival (S) 
ms.ddl$S$AGE<- ifelse(ms.ddl$S$Age > 0, "Adult", "Juvenile")
#head(ms.ddl$S)

#so young are classified as ageclass=1 and adults as ageclass=0
ms.ddl$S$ageclass=0
ms.ddl$S$ageclass [ms.ddl$S$age==0]=1

#make dummy variables for juv/adult (we didn´t use it)
#ms.ddl$S$adult=0
#ms.ddl$S$juv=0
#ms.ddl$S$adult [ms.ddl$S$ageclass==0]=1
#ms.ddl$S$juv [ms.ddl$S$ageclass==1]=1

# Add age classes (juvenile, adult) for transition (Psi) parameters
ms.ddl$Psi$AGE<- ifelse(ms.ddl$Psi$Age > 0, "Adult", "Juvenile")
#head(ms.ddl$Psi)



# Add invasion as a covariate for survival (S) 
inv <- read.csv("tsi_by_region_for_ss_96_23_50perc.csv")
inv$invaded <- ifelse(inv$TSI < 0, 0, 1)
head(inv)
colnames(inv) <- c("ref", "time", "Pop", "Year_invaded", "TSI", "stratum", "invaded")
inv$time <- as.factor(inv$time)
inv$stratum <- as.factor(inv$stratum)
inv$invaded <- as.factor(inv$invaded)


#Add SPEI
gridmet.spei.inv <- read.csv("SPEI_per_yr_long_extreme_inv.csv")

gridmet.spei.inv$year <- as.factor(gridmet.spei.inv$year)

inv.spei <- inv %>%
  left_join(gridmet.spei.inv, by = c("time" = "year", "Pop" = "Pop")) %>%
  mutate(n_spei = case_when(spei == "Normal" ~ 1,
                            spei == "Drought" ~ 2,
                            spei == "Wet" ~ 3)) %>% 
  dplyr::select(-X,
                -TSI,
                -Year_invaded)

head(inv.spei)

inv.spei$spei <- as.factor(inv.spei$spei)

ms.ddl$S <- inner_join(ms.ddl$S, inv.spei)
head(ms.ddl$S)
ms.ddl$p <- inner_join(ms.ddl$p, inv.spei)
head(ms.ddl$p)
ms.ddl$Psi <- inner_join(ms.ddl$Psi, inv.spei)
head(ms.ddl$Psi)


#Add a dummy variable for droughts
ms.ddl$S$drought = 0
ms.ddl$S$drought [ms.ddl$S$n_spei == 2] = 1


# Convert ms.ddl$S$AGE, ms.ddl$S$spei and ms.ddl$S$drought to factor with appropriate levels
ms.ddl$S$AGE <- factor(ms.ddl$S$AGE, levels = c("Juvenile", "Adult"))
ms.ddl$S$spei <- factor(ms.ddl$S$spei, levels = c("Normal", "Drought", "Wet"))

#ms.ddl$S$adult <- factor(ms.ddl$S$adult)
#ms.ddl$S$juv <- factor(ms.ddl$S$juv)

ms.ddl$S$drought = factor(ms.ddl$S$drought)

str(ms.ddl$S)

# Add wetland type (lacustrine vs palustrine) as a detection parameter
ms.ddl$p$wt <- ifelse(ms.ddl$p$stratum %in% c(1:3),"L", "P")
head(ms.ddl$p)

#### Add distance (m) between populations centroids as a covariate for Psi 
#From Jeff Laake: http://www.phidot.org/forum/viewtopic.php?f=21&t=643

ms.ddl$Psi$distance=0.000001 #value close to cero so we can include log(distance) into the model
ms.ddl$Psi$distance[ms.ddl$Psi$stratum=="1"&ms.ddl$Psi$tostratum=="2"]=214911
ms.ddl$Psi$distance[ms.ddl$Psi$stratum=="1"&ms.ddl$Psi$tostratum=="3"]=326151
ms.ddl$Psi$distance[ms.ddl$Psi$stratum=="1"&ms.ddl$Psi$tostratum=="4"]=251481
ms.ddl$Psi$distance[ms.ddl$Psi$stratum=="1"&ms.ddl$Psi$tostratum=="5"]=449329
ms.ddl$Psi$distance[ms.ddl$Psi$stratum=="1"&ms.ddl$Psi$tostratum=="6"]=357643

ms.ddl$Psi$distance[ms.ddl$Psi$stratum=="2"&ms.ddl$Psi$tostratum=="1"]=214911
ms.ddl$Psi$distance[ms.ddl$Psi$stratum=="2"&ms.ddl$Psi$tostratum=="3"]=111314
ms.ddl$Psi$distance[ms.ddl$Psi$stratum=="2"&ms.ddl$Psi$tostratum=="4"]=57937
ms.ddl$Psi$distance[ms.ddl$Psi$stratum=="2"&ms.ddl$Psi$tostratum=="5"]=238638
ms.ddl$Psi$distance[ms.ddl$Psi$stratum=="2"&ms.ddl$Psi$tostratum=="6"]=147165

ms.ddl$Psi$distance[ms.ddl$Psi$stratum=="3"&ms.ddl$Psi$tostratum=="1"]=326151
ms.ddl$Psi$distance[ms.ddl$Psi$stratum=="3"&ms.ddl$Psi$tostratum=="2"]=111314
ms.ddl$Psi$distance[ms.ddl$Psi$stratum=="3"&ms.ddl$Psi$tostratum=="4"]=95620
ms.ddl$Psi$distance[ms.ddl$Psi$stratum=="3"&ms.ddl$Psi$tostratum=="5"]=132636
ms.ddl$Psi$distance[ms.ddl$Psi$stratum=="3"&ms.ddl$Psi$tostratum=="6"]=58575

ms.ddl$Psi$distance[ms.ddl$Psi$stratum=="4"&ms.ddl$Psi$tostratum=="1"]=251481
ms.ddl$Psi$distance[ms.ddl$Psi$stratum=="4"&ms.ddl$Psi$tostratum=="2"]=57937
ms.ddl$Psi$distance[ms.ddl$Psi$stratum=="4"&ms.ddl$Psi$tostratum=="3"]=95620
ms.ddl$Psi$distance[ms.ddl$Psi$stratum=="4"&ms.ddl$Psi$tostratum=="5"]=227804
ms.ddl$Psi$distance[ms.ddl$Psi$stratum=="4"&ms.ddl$Psi$tostratum=="6"]=107894

ms.ddl$Psi$distance[ms.ddl$Psi$stratum=="5"&ms.ddl$Psi$tostratum=="1"]=449329
ms.ddl$Psi$distance[ms.ddl$Psi$stratum=="5"&ms.ddl$Psi$tostratum=="2"]=238638
ms.ddl$Psi$distance[ms.ddl$Psi$stratum=="5"&ms.ddl$Psi$tostratum=="3"]=132636
ms.ddl$Psi$distance[ms.ddl$Psi$stratum=="5"&ms.ddl$Psi$tostratum=="4"]=227804
ms.ddl$Psi$distance[ms.ddl$Psi$stratum=="5"&ms.ddl$Psi$tostratum=="6"]=142605


ms.ddl$Psi$distance[ms.ddl$Psi$stratum=="6"&ms.ddl$Psi$tostratum=="1"]=357643
ms.ddl$Psi$distance[ms.ddl$Psi$stratum=="6"&ms.ddl$Psi$tostratum=="2"]=147165
ms.ddl$Psi$distance[ms.ddl$Psi$stratum=="6"&ms.ddl$Psi$tostratum=="3"]=58575
ms.ddl$Psi$distance[ms.ddl$Psi$stratum=="6"&ms.ddl$Psi$tostratum=="4"]=107894
ms.ddl$Psi$distance[ms.ddl$Psi$stratum=="6"&ms.ddl$Psi$tostratum=="5"]=142605

#Original code from Meghan used the scaled distance but we won´t do that because we include log(distance) into the model
#ms.ddl$Psi$distance <- scale(ms.ddl$Psi$distance)

# Fix survival to 0 for stratum 1 (Paynes Prairie) for all years before 2018 (no snail kites)
ms.ddl$S$fix<- ifelse((ms.ddl$S$stratum == '1' & !(ms.ddl$S$time %in% c('2018','2019','2020','2021', '2022', '2023'))),
                      0,
                      NA)
# Fix dispersal to 0 for stratum 1 (Paynes Prairie) for all years before 2018 (no snail kites)
ms.ddl$Psi$fix<- ifelse((ms.ddl$Psi$stratum == '1' & !(ms.ddl$Psi$time %in% c('2018','2019','2020','2021', '2022', '2023')) |
                           ms.ddl$Psi$tostratum == '1' & !(ms.ddl$Psi$time %in% c('2018','2019','2020','2021', '2022', '2023'))),
                        0,
                        NA)
# Fix p to 0 for stratum 1 (Paynes Prairie) for all years before 2018 (no snail kites)
ms.ddl$p$fix<- ifelse(ms.ddl$p$stratum == '1' & !(ms.ddl$p$time %in% c('2018','2019','2020','2021', '2022', '2023')),
                      0,
                      NA)


## Prepare multistate model #####
run.ms=function() {
  
  #detection by wetland type (lacustrine vs palustrine)
  p.wetland = list(formula = ~wt) #detection varies by wetland type
  
  #adding time and stratum as co-variates resulted in better models (based on AIC) but had convergence issues:
  
  #original from Josh (2020 PVA): 
  #p.time.stratum = list(formula = ~time + stratum)
  
  #p.stratum = list(formula = ~  stratum)
  
  
  # state transitions 
  #adding both log(distance) and AGE resulted in a better model than only distance or only AGE. Adding stratum:tostratum had convergence issues.
  
  #Psi.dist.age =	list(formula = ~distance +  AGE)  
  Psi.logDist.age =	list(formula = ~log(distance) +  AGE)
  #Psi.logDist =	list(formula = ~log(distance) )
  #Psi.age =	list(formula = ~  AGE)
  
  #original from Josh (2020 PVA):
  #Psi.stratum  =	list(formula = ~-1 + stratum:tostratum)
  
  # survival
  
  #speiXPop interaction: model with large SE due to lack of convergence
  #S.strat.age.speiXpop.inv  =	list(formula = ~AGE + stratum *spei + invaded) 
  
  #What Josh did for 2020 PVA (results were unexpected - with higher survival in droughts)
  #S.stratum.inv.speiXage = list(formula = ~AGE + stratum + spei + AGE:spei+invaded)
  
  #trials that didn´t converge:
  #S.stratum.inv.speiXAge.juvXPop = list(formula = ~ stratum + spei + adult:spei+ juv:spei:stratum + invaded)
  #S.stratum.inv.speiXAge.juvXPop2 = list(formula = ~ stratum + spei + AGE + AGE:spei + juv:stratum + invaded)
  #S.stratum.inv.spei.Age.juvXPop = list(formula = ~ stratum + spei + adult + juv:stratum + invaded)
  #S.stratum.inv.spei.Age.juvXPop2 = list(formula = ~  spei + adult + juv:stratum + invaded)
  #S.stratum.inv.droughtXPop.droughtXAge = list(formula = ~ stratum + AGE + AGE:drought + drought:stratum + invaded)
  
  S.stratum.inv.droughtXPop.Age = list(formula = ~ stratum + AGE + drought:stratum + invaded)
  
  #simpler models
  #S.stratum.droughtXPop.Age2 = list(formula = ~ stratum + AGE + drought*stratum)
  #S.strat.age.spei.inv.pop  =	list(formula = ~AGE + stratum  + invaded + spei)
  
  ms.model.list=create.model.list("Multistrata")
  ms.results=mark.wrapper(ms.model.list,
                          data=ms.pr, ddl=ms.ddl, adjust=FALSE, output=FALSE)
  return(ms.results)
}

#Run multistate model, takes ~30 minutes

system.time(ms.results_inv.droughtXPop.Age<- run.ms())

#save.image("ms.results_inv.droughtXPop.Age.RData")

# predictions for survival -------------------------------------------------------------
load("ms.results_inv.droughtXPop.Age.RData")

# Beta values

beta_df <- ms.results_inv.droughtXPop.Age[[1]]$results$beta

# Extract beta estimates and set names
beta <- beta_df$estimate
names(beta) <- rownames(beta_df)

# Extract variance-covariance matrix
vcv <- ms.results_inv.droughtXPop.Age[[1]]$results$beta.vcv

# Extract factor levels from the design data to ensure consistency
AGE_levels <- unique(ms.ddl$S$AGE)
stratum_levels <- unique(ms.ddl$S$stratum)
invaded_levels <- unique(ms.ddl$S$invaded)
#spei_levels <- unique(ms.ddl$S$spei)
drought_levels <- unique(ms.ddl$S$drought)

# Create combinations for prediction
pred_data <- expand.grid(
  AGE = AGE_levels,
  stratum = stratum_levels,
 invaded = invaded_levels,
#spei_cat = spei_levels
  drought = drought_levels
)

# Ensure variables are factors with correct levels
pred_data$AGE <- factor(pred_data$AGE, levels = levels(ms.ddl$S$AGE))
pred_data$stratum <- factor(pred_data$stratum, levels = levels(ms.ddl$S$stratum))
pred_data$invaded <- factor(pred_data$invaded, levels = levels(ms.ddl$S$invaded))
#pred_data$spei <- factor(pred_data$spei, levels = levels(ms.ddl$S$spei))
pred_data$drought <- factor(pred_data$drought, levels = levels(ms.ddl$S$drought))

# Generate the design matrix using the model formula
design <- model.matrix(~ stratum + AGE + drought:stratum + invaded, data = pred_data)

#check consistency between design and beta values

# Get indices of 'S' parameters
indices_S <- grep("^S:", names(beta))
# Extract the beta coefficients for 'S'
beta_S <- beta[indices_S]
# Verify length
length(beta_S)  
# Number of columns in the design matrix
ncol(design)
# Names of beta coefficients
names(beta_S)
# Column names of the design matrix
colnames(design)

#check if it is needed to change the reference values to match colnames(beta_S)
# Relevel AGE to have 'Adult' as the reference level
#pred_data$AGE <- relevel(pred_data$AGE, ref = "adult")


# Adjust design matrix column names to match beta_S names
colnames(design) <- paste0("S:", colnames(design))

#there is no data for this combination (PP and drought)
design <- design[, -which(colnames(design) == "S:stratum1:drought1")]

# Check that column names match
all(colnames(design) == names(beta_S))
# Should return TRUE

# Extract the variance-covariance matrix for 'S' parameters
vcv_S <- vcv[indices_S, indices_S]


# Calculate eta (linear predictor)
eta <- design %*% beta_S

# Compute variance of eta
var_eta <- diag(design %*% vcv_S %*% t(design))
se_eta <- sqrt(var_eta)

# Transform to survival probabilities
alpha <- 0.05
z <- qnorm(1 - alpha / 2)
eta_lcl <- eta - z * se_eta
eta_ucl <- eta + z * se_eta

S <- plogis(eta)
S_lcl <- plogis(eta_lcl)

S_ucl <- plogis(eta_ucl)


# Assemble results
pred_data$estimate <- as.vector(S)
pred_data$lcl <- as.vector(S_lcl)
pred_data$ucl <- as.vector(S_ucl)
pred_data$se <- (pred_data$estimate- pred_data$lcl)/z

#change stratum numbers to names
stratum_labels <- c("1" = "PP", "2" = "KRV", "3" = "OKEE", "4" = "SJM", "5" = "EVER", "6" = "EAST")

pred_data$stratum <- factor(pred_data$stratum,
                            levels = names(stratum_labels),
                            labels = stratum_labels)

#change invaded numbers to names
pred_data$invaded <- factor(pred_data$invaded, 
                            levels = c(0, 1), 
                            labels = c( "No","Yes"))

pred_data$drought <- factor(pred_data$drought, 
                            levels = c(0, 1), 
                            labels = c( "No","Yes"))


write.csv(pred_data, "data/survival_estimates_PVA_age.droughtXpop.inv.csv")


### plots -------------------------------------------------------------------

#change the order of categories for plotting
pred_data$stratum <- factor(pred_data$stratum, levels = c("EAST", "EVER", "KRV", "OKEE", "PP", "SJM"))

#pred_data$spei <- factor(pred_data$spei, 
#                         levels = c("Drought","Normal", "Wet"))

# Plot survival estimates
ggplot(pred_data, aes(x = drought, y = estimate, color = invaded)) +
  geom_point(position = position_dodge(width = 0.2), size =3) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.2, position = position_dodge(width = 0.2)) +
  facet_grid(AGE ~ stratum) +
  labs(
    # title = "Predicted Survival Probabilities",
    x = "Droughts",
    y = "Survival",
    color = "Invaded"
  ) +
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))


#Both ages in the same plot
ggplot(pred_data, aes(x = drought, color = AGE, shape=invaded)) +
  geom_point(aes(y = estimate), size = 1.5) +
  geom_linerange(aes(ymin = lcl, ymax = ucl), alpha = 0.5) +
  scale_color_manual("", values = wes_palette("Cavalcanti1", 2)) +
  #coord_flip() +
  theme_bw() +
  facet_wrap(~ stratum)

#final plot (only for invaded = Yes)
pred_data_inv <- pred_data %>% 
  filter(invaded == "Yes")

ggplot(pred_data_inv, aes(x = drought, y = estimate, color= AGE)) +
  geom_point(position = position_dodge(width = 0.2), size =3) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.2, position = position_dodge(width = 0.2)) +
  facet_grid( ~ stratum) +
  labs(
    # title = "Predicted Survival Probabilities",
    x = "Droughts",
    y = "Apparent survival (CI)",
    color = "Age"
  ) +
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))


# predictions for transition ----------------------------------------------

top1<-ms.results_inv.droughtXPop.Age[[1]]
design=top1$design.matrix

top1$results$beta

xx=get.real(top1,"Psi",vcv=TRUE) 
#takes ~20 minutes

#Transition matrix will be estimated per year but all years will have the same values prior 2018 and after 2018 (when PP started having Snail kite occurrences). I decided to export the transition matrix for any year after 2018 given that´s the current scenario.

#Adults matrix (has to be after 1997)
time=2022
res_ad=TransitionMatrix(xx$estimates[xx$estimates$time==time&xx$estimates$age==1,],vcv=xx$vcv.real)

trans_matrix_ad <- res_ad$TransitionMat
colnames(trans_matrix_ad) <- c("PP", "KRV", "OKEE", "SJM", "EVER", "EAST")
rownames(trans_matrix_ad) <- c("PP", "KRV", "OKEE", "SJM", "EVER", "EAST")

se.trans_matrix_ad <- res_ad$se.TransitionMat
colnames(se.trans_matrix_ad) <- c("PP", "KRV", "OKEE", "SJM", "EVER", "EAST")
rownames(se.trans_matrix_ad) <- c("PP", "KRV", "OKEE", "SJM", "EVER", "EAST")

#juveniles matrix
time=2022
res_juv=TransitionMatrix(xx$estimates[xx$estimates$time==time&xx$estimates$age==0,],vcv=xx$vcv.real)

trans_matrix_juv <- res_juv$TransitionMat
colnames(trans_matrix_juv) <- c("PP", "KRV", "OKEE", "SJM", "EVER", "EAST")
rownames(trans_matrix_juv) <- c("PP", "KRV", "OKEE", "SJM", "EVER", "EAST")

se.trans_matrix_juv <- res_juv$se.TransitionMat
colnames(se.trans_matrix_juv) <- c("PP", "KRV", "OKEE", "SJM", "EVER", "EAST")
rownames(se.trans_matrix_juv) <- c("PP", "KRV", "OKEE", "SJM", "EVER", "EAST")

write.csv(trans_matrix_ad, "data/trans_matrix_adults.csv")
write.csv(trans_matrix_juv, "data/trans_matrix_juveniles.csv")
write.csv(se.trans_matrix_ad, "data/se.trans_matrix_adults.csv")
write.csv(se.trans_matrix_juv, "data/se.trans_matrix_juveniles.csv")




