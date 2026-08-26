library(ggplot2)
library(dplyr)
library(gratia)
library(glmmTMB)
library(MuMIn)

# ESTIMATING NEST SUCCESS PROBABILITY AND NUMBER OF YOUNG FLEDGED #### 

#see README for column descriptions

# NEST SUCCESS #### 

#load data
nests <- read.csv("data/nest_survival_data_Feb2024.csv")
head(nests)

#gridmet.spei <- read.csv("data/SPEI_per_yr_long.csv")
#gridmet.spei <- read.csv("data/spei_per_year_avg02-07_long.csv")
gridmet.spei.inv <- read.csv("data/SPEI_per_yr_long_extreme_inv.csv")

#check which combinations do not exist
all_comb <- expand.grid(spei = c("Drought", "Normal", "Wet"), Pop = c("EAST", "EVER", "OKEE", "PP", "SJM", "KRV"), inv.cat = c("Yes", "No"))


gridmet.spei_inv_NA <- gridmet.spei.inv %>% 
  right_join(all_comb) %>% 
  filter(is.na(year))

#EAST-Drought-Yes
#PP-Drought-Yes
#SJM-Drought-Yes
#PP-Wet-No

nest.spei<- nests %>%
  left_join(gridmet.spei.inv, join_by(Year==year, Pop == Pop))

#I filtered PP-Wet because the nests haven´t failed yet and the prediction was 1 for nest success
nest.spei_filt <- nest.spei %>% 
  filter(Pop != "PP"  | spei != "Wet")

nest.spei_filt$spei <- as.factor(nest.spei_filt$spei)
nest.spei_filt$inv.cat <- as.factor(nest.spei_filt$inv.cat)
nest.spei_filt$Pop <- as.factor(nest.spei_filt$Pop)

#let´s check which combinations have no nests
nest.spei_summary <- nest.spei %>% 
  group_by(spei, Pop, inv.cat) %>% 
  summarise(count = n())

comb_no_nests <- all_comb %>% 
  left_join(nest.spei_summary) %>% 
  filter(is.na(count))



## models ------------------------------------------------------------------
# I want "Normal" to be in the intercept
nest.spei_filt$spei <- relevel(nest.spei_filt$spei, ref = "Normal")

#triple interaction
#ns.speiXinvXpop <- glmmTMB(FateInt ~   spei*inv.cat*Pop, offset = log(exposure),
#                           family=binomial(link="cloglog"), data=nest.spei_filt)

#summary(ns.speiXinvXpop)

#additive effect of "invaded" to simplify (because there are many combinations with NA in the previous model)
ns.speiXpop.inv <- glmmTMB(FateInt ~   Pop*spei+inv.cat, offset = log(exposure),
                            family=binomial(link="cloglog"), data=nest.spei_filt)

summary(ns.speiXpop.inv)

#with no interaction (because there are combinations with NA in the previous model)
#ns.spei.pop.inv <- glmmTMB(FateInt ~   Pop+spei+inv.cat, offset = log(exposure),
#                           family=binomial(link="cloglog"), data=nest.spei_filt)


#summary(ns.spei.pop.inv)

#just to check if our models are better than the null model
ns.null <- glmmTMB(FateInt ~   1, offset = log(exposure),
                           family=binomial(link="cloglog"), data=nest.spei_filt)


model_list <- list(ns.spei.pop.inv = ns.spei.pop.inv, ns.speiXpop.inv = ns.speiXpop.inv, ns.speiXinvXpop = ns.speiXinvXpop,  ns.null = ns.null)


model_sel <- model.sel(model_list, rank = "AICc")



### predictions -------------------------------------------------------------

#model ns.speiXpop.inv

# Create ns.newdata with all combinations
ns.newdata <- expand.grid(
  spei = levels(nest.spei$spei),
  Pop = levels(nest.spei$Pop),
  inv.cat = levels(nest.spei$inv.cat)
)

ns.newdata$exposure <- 1 #set exposure days to 1 for prediction
dsr_pred <- predict(ns.speiXpop.inv, newdata = ns.newdata, re.form = NULL, type="response", se.fit = T)
ns.newdata$pred <- dsr_pred$fit
ns.newdata$DSR <- 1- dsr_pred$fit
ns.newdata$DSR_SE <- dsr_pred$se.fit
ns.newdata$DSR_lower <- (ns.newdata$DSR - (1.96 * dsr_pred$se.fit))
ns.newdata$DSR_upper <- (ns.newdata$DSR + (1.96 * dsr_pred$se.fit))
ns.newdata$ns <- (ns.newdata$DSR)^58 #convert to nest success by raising to length of nest cycle (58 days)
ns.newdata$ns_upper <- ns.newdata$DSR_upper^58
ns.newdata$ns_lower <- ns.newdata$DSR_lower^58
ns.newdata$NS_SE <- (ns.newdata$ns- ns.newdata$ns_lower)/1.96


# replace values to NA for combinations with no nests.
#PP had no droughts since snail kites are there so we kept the estimated nest success by the model
nest.spei_summary2 <- nest.spei %>% 
  group_by(spei, Pop) %>% 
  summarise(count = n())

comb_no_nests <- all_comb %>% 
  left_join(nest.spei_summary2) %>% 
  filter(is.na(count))


ns.newdata_NA <- ns.newdata %>%
  mutate(
    pred = ifelse(Pop %in% c( "SJM", "EAST") & spei == "Drought", NA, pred),
    DSR = ifelse(Pop %in% c( "SJM", "EAST") & spei == "Drought", NA, DSR),
    DSR_SE = ifelse(Pop %in% c( "SJM", "EAST") & spei == "Drought", NA, DSR_SE),
    DSR_lower = ifelse(Pop %in% c( "SJM", "EAST") & spei == "Drought", NA, DSR_lower),
    DSR_upper = ifelse(Pop %in% c( "SJM", "EAST") & spei == "Drought", NA, DSR_upper),
    ns = ifelse(Pop %in% c( "SJM", "EAST") & spei == "Drought", NA, ns),
    ns_lower = ifelse(Pop %in% c( "SJM", "EAST") & spei == "Drought", NA, ns_lower),
    ns_upper = ifelse(Pop %in% c( "SJM", "EAST") & spei == "Drought", NA, ns_upper),
    NS_SE = ifelse(Pop %in% c( "SJM", "EAST") & spei == "Drought", NA, NS_SE)
  )

ns.newdata_NA2 <- ns.newdata_NA %>% 
  left_join(nest.spei_summary) %>% 
  dplyr::rename(Nnest = count)
  
write.csv(ns.newdata_NA2, "data/PVA_nest_success_estimates.csv")

#plot

ns.newdata_NA$spei <- factor(ns.newdata_NA$spei, 
                            levels = c( "Drought","Normal", "Wet"))

ns.newdata_NA %>% 
  ggplot(aes(x = spei, y = ns, color= inv.cat)) +
  geom_errorbar( aes(x = spei, ymin= ns_lower, ymax= ns_upper), width = 0.2, 
                 position = position_dodge(width = 0.5))+
  geom_point(data= ns.newdata_NA,position = position_dodge(width = 0.5),size=3)+
  facet_grid(~ Pop) +
  ylim(c(0,1)) +
  labs(title = '',
       y = "Nest success probability", 
       x = 'SPEI Category',
       color = "Invaded") +
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

#final plot (only invaded scenario)

ns.newdata_NA %>% 
  filter(inv.cat == "Yes") %>% 
  ggplot(aes(x = spei, y = ns)) +
  geom_errorbar( aes(x = spei, ymin= ns_lower, ymax= ns_upper), width = 0.2)+
  geom_point(size=3)+
  facet_grid(~ Pop) +
  ylim(c(0,1)) +
  labs(title = '',
       y = "Nest success probability (CI)", 
       x = 'SPEI Category') +
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))



### residuals ---------------------------------------------------------------

# Pearson residuals
pearson_residuals <- residuals(ns.speiXpop.inv, type = "pearson")

# Deviance residuals
deviance_residuals <- residuals(ns.speiXpop.inv, type = "deviance")

# Response residuals
response_residuals <- residuals(ns.speiXpop.inv, type = "response")

# Fitted values
fitted_values <- predict(ns.speiXpop.inv, type = "response")

# Plot Pearson residuals
plot(fitted_values, pearson_residuals, 
     ylab = "Pearson Residuals", xlab = "Fitted values",
     main = "Residuals vs Fitted values")
abline(h = 0, col = "red")

# Histogram of residuals
hist(pearson_residuals, main = "Histogram of Pearson Residuals", xlab = "Residuals")

#checking for overdispersion
# Sum of squared Pearson residuals
sum_pearson <- sum(pearson_residuals^2)

# Residual degrees of freedom
df_residual <- df.residual(ns.speiXpop.inv)

# Overdispersion ratio
overdispersion_ratio <- sum_pearson / df_residual

overdispersion_ratio #good

# DHARMa package for uniform residuals
library(DHARMa)

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = ns.speiXpop.inv)

# Plot diagnostics from DHARMa
plot(sim_res)

# YOUNG FLEDGED ####

prod <- read.csv("data/yng_fledged_data_Feb2024.csv")
head(prod)

# Check the distribution of the response variable
hist(prod.spei$Fledglings, main = "Histogram of Fledglings", xlab = "Fledglings")

prod.spei<- prod  %>%
  left_join(gridmet.spei.inv, join_by(Year==year, Pop == Pop)) 

str(prod.spei)
prod.spei$spei <- as.factor(prod.spei$spei)
prod.spei$inv.cat <- as.factor(prod.spei$inv.cat)
prod.spei$Pop <- as.factor(prod.spei$Pop)

prod_summary <- prod.spei %>% 
  group_by(spei, Pop) %>% 
  summarise(count_nests=n())

#No nests:
#EAST-Drought
#PP-Drought
#SJM-Drought

## models ------------------------------------------------------------------
prod.spei$spei <- relevel(prod.spei$spei, ref = "Normal")

#p.speiXinvXpop <- glmmTMB(Fledglings ~ Pop*spei*inv.cat, family = poisson, data = prod.spei)

p.speiXpop.inv <- glmmTMB(Fledglings ~ Pop*spei+inv.cat, family = poisson, data = prod.spei)

summary(p.speiXpop.inv)

p.spei.inv.pop <- glmmTMB(Fledglings ~ Pop+spei+inv.cat, family = poisson, data = prod.spei)

summary(p.spei.inv.pop)

p.null <- glmmTMB(Fledglings ~ 1, family = poisson, data = prod.spei)


model_list.p <- list(p.spei.inv.pop = p.spei.inv.pop, p.speiXpop.inv = p.speiXpop.inv, p.speiXinvXpop = p.speiXinvXpop, p.null = p.null)

model_sel.p <- model.sel(model_list.p, rank = "AICc")
#the model with no interactions is a better model. 

### predictions -------------------------------------------------------------

# Create ns.newdata with all combinations
p.newdata <- expand.grid(
  spei = levels(prod.spei$spei),
  Pop = levels(prod.spei$Pop),
  inv.cat = levels(prod.spei$inv.cat)
)


prod_pred <- predict(p.spei.inv.pop, newdata = p.newdata, type="response", se.fit = T, re.form = NA)
p.newdata$Yng_cnt <- prod_pred$fit
p.newdata$lower <- (p.newdata$Yng_cnt - (1.96 * prod_pred$se.fit))
p.newdata$upper <- (p.newdata$Yng_cnt + (1.96 * prod_pred$se.fit))
p.newdata$se <- prod_pred$se.fit

# replace values to NA for combiantions with no nests
p.newdata_NA <- p.newdata %>%
  mutate(
    Yng_cnt = ifelse(Pop %in% c( "SJM", "EAST") & spei == "Drought", NA, Yng_cnt),
    lower = ifelse(Pop %in% c( "SJM", "EAST") & spei == "Drought", NA, lower),
    upper = ifelse(Pop %in% c( "SJM", "EAST") & spei == "Drought", NA, upper),
    se = ifelse(Pop %in% c( "SJM", "EAST") & spei == "Drought", NA, se)
  )

Nnests <- prod.spei %>% 
  group_by(spei, Pop, inv.cat) %>% 
  summarise(Nnests=n())

p.newdata_Nnest <- p.newdata_NA %>% 
  left_join(Nnests)

write.csv(p.newdata_Nnest, "data/PVA_young_fledged_estimates.csv")

# plot

p.newdata_NA$spei <- factor(p.newdata_NA$spei, 
                          levels = c( "Drought","Normal", "Wet"))

p.newdata_NA %>% 
  filter(inv.cat == "Yes") %>% 
  ggplot(aes(x = spei, y = Yng_cnt)) +
  geom_errorbar( aes(x = spei, ymin= lower, ymax= upper), width=0.2)+
  geom_point(size=3)+
  facet_grid(~ Pop) +
  labs(title = '',
       y = "Number of young fledged (CI)", 
       x = 'SPEI Category')+
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))


### residuals ---------------------------------------------------------------

# Extract residuals
pearson_residuals <- residuals(p.spei.inv.pop, type = "pearson")
deviance_residuals <- residuals(p.spei.inv.pop, type = "deviance")
response_residuals <- residuals(p.spei.inv.pop, type = "response")

# Extract fitted values
fitted_values <- predict(p.spei.inv.pop, type = "response")

# Plot Pearson residuals
plot(fitted_values, pearson_residuals, 
     ylab = "Pearson Residuals", xlab = "Fitted values",
     main = "Pearson Residuals vs Fitted Values")
abline(h = 0, col = "red")

# Plot Deviance residuals
plot(fitted_values, deviance_residuals, 
     ylab = "Deviance Residuals", xlab = "Fitted values",
     main = "Deviance Residuals vs Fitted Values")
abline(h = 0, col = "blue")

# Histogram of Pearson Residuals
hist(pearson_residuals, main = "Histogram of Pearson Residuals", xlab = "Pearson Residuals")

# Histogram of Deviance Residuals
hist(deviance_residuals, main = "Histogram of Deviance Residuals", xlab = "Deviance Residuals")

# Check for overdispersion
sum_pearson <- sum(pearson_residuals^2)
df_residual <- df.residual(p.spei.inv.pop)
overdispersion_ratio <- sum_pearson / df_residual
print(overdispersion_ratio)

# Use DHARMa for diagnostic plots
sim_res <- simulateResiduals(fittedModel = p.spei.inv.pop)
plot(sim_res)
