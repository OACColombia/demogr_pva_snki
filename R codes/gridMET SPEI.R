
### Download gridMET observed climate data for SPEI, SPI, and PDSI ###
## For comparison against MACA projections of SPEI ##

remotes::install_github("mikejohnson51/AOI") # suggested!
remotes::install_github("mikejohnson51/climateR")

library(tidyverse)
library(dplyr)
library(sf)                 # dealing with vector data
library(exactextractr)      # Fast zonal statistics
library(climateR)           # Download gridMET observational climate data
library(rnaturalearth)      # Access state shapefiles
library(rnaturalearthdata)  # Access state shapefiles
library(rnaturalearthhires) # Access state shapefiles
library(SPEI)               # Calculate SPEI
library(furrr)              # Parallel processing
library(future)             # Parallel processing
library(tictoc)             # Measure time for code to run
library(wesanderson)        # Color palettes
library(terra)
library(AOI)



source('R codes/functions updated.R')


# Load MU shp
mu_Meghan <- st_read ("./data", "MUs_Meghan")

# Access shapefiles of states
states<- ne_states(country = "United States of America", returnclass = "sf") %>% 
  st_transform(crs = 32617)

# Define color palette for MUs
pal<- c(wes_palette('Darjeeling1', n = 5),
        wes_palette('Darjeeling2', n = 2)[2])
pal<- pal[c(1,4,3,2,5,6)]  #reorder

bounds<- st_sfc(st_point(c(-83,24.5)),
                st_point(c(-80,30.5)),
                crs = 4326) %>% 
  st_transform(crs = 32617)

ggplot() +
  geom_sf(data = states) +
  geom_sf(data = mu_Meghan, aes(fill = name)) +
  scale_fill_manual("", values = pal) +
  coord_sf(xlim = st_bbox(bounds)[c("xmin","xmax")],
           ylim = st_bbox(bounds)[c("ymin","ymax")],
           crs = 32617) +
  theme_bw() +
  theme(panel.grid = element_blank())



# Dissolve boundaries across neighboring sites within same MU

mu <- mu_Meghan %>% 
  group_by(name) %>%
  summarize() %>%  #this is the function that performs a spatial union
  st_transform(crs = 4326)  #transform to match CRS used in getMACA()

#Download gridMET data (4 km res; daily res)
mu.list<- split(mu, seq(nrow(mu)))   #split each MU into separate elements


# Jan 1994 - Dec 2025 (Observed climate data)
tic()
gridmet.res <- map(mu.list,
                  ~gridMET_summary(AOI = .x,
                                    varnames = c("pr",'tmmn','tmmx'),
                                    startDate = "1994-01-01", #start 1994 for 18 month SPEI
                                    endDate = "2025-12-31")
)
toc()


## Label MUs for each dataset and merge together
names(gridmet.res)<- mu$name
#names(gridmet.res.pdsi)<- mu2$MU_Name

#convert from list to DF
gridmet.res<- gridmet.res %>% 
  bind_rows(.id = "MU")



# calculate midpoint of temp and convert from K -> C for SPEI
gridmet<- gridmet.res %>% 
  mutate(tmid = (tmax + tmin)/2,
         tmidC = tmid - 273.15,
         .before = tmax)

# export data for ease of access
# write.csv(gridmet, "data/SNKI MU 1994-2025 gridMET data.csv", row.names = F)
# gridmet<- read.csv("data/SNKI MU 1994-2025 gridMET data.csv", as.is = T)

# add year-date for plotting
gridmet$my<- lubridate::as_date(paste(gridmet$year, gridmet$month, 01, sep = "-"))


# Plot precip
ggplot(gridmet, aes(my, prcp, color = MU)) +
  geom_line(size = 1) +
  scale_color_manual("", values = pal) +
  labs(x = "Year", y = "Monthly Precipitation (mm)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14, face = "bold"),
        strip.background = element_blank(),
        legend.text = element_text(size = 12))

# Plot temp
ggplot(gridmet, aes(my, tmidC, color = MU)) +
  geom_line(size = 1) +
  scale_color_manual("", values = pal) +
  labs(x = "Year", y = "Monthly Avg Temperature (\u00B0C)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14, face = "bold"),
        strip.background = element_blank(),
        legend.text = element_text(size = 12))



## Calculate SPEI (Standardized Precipitation Evapotranspiration Index)

#calculate potential evapotranspiration based on thornthwaite method
avg.lat<- map(mu.list, function(x) mean(st_bbox(x)[c("ymin","ymax")])) %>% 
  set_names(mu$name)

gridmet2<- gridmet %>% 
  split(.$MU) %>%
  map2(.x = ., .y = avg.lat, .f = ~{.x %>% 
      mutate(PET = as.vector(thornthwaite(tmidC, .y)))
  }) %>% 
  bind_rows()

ggplot(gridmet2, aes(my, PET, color = MU)) +
  geom_line(size = 1) +
  scale_color_manual("", values = pal) +
  labs(x = "Year", y = "Monthly PET (mm)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14, face = "bold"),
        strip.background = element_blank(),
        legend.text = element_text(size = 12))



# Calculate balance (avg precip - potential evapotranspiration) as input for SPEI
gridmet2<- gridmet2 %>% 
  mutate(BAL = prcp - PET)


#calculate 18-month SPEI with log-Logistic fitting function using equal weights across months
gridmet.spei<- gridmet2 %>% 
  split(.$MU) %>% 
  map(., ~spei(ts(.[,'BAL'], freq = 12, start = 1994), scale = 18, na.rm = T))



# Reformat SPEI data so it is more user friendly and define SPEI categories

#7 categories:
#gridmet.spei.df<- gridmet.spei %>% 
#  map2(.x = ., .y = split(gridmet2, gridmet2$MU), ~{
#    .$fitted %>% 
#      data.frame() %>% 
#      rename(index = ".") %>% 
#      mutate(year = .y$year, month = .y$month)
#  }) %>% 
#  bind_rows(.id = "MU") %>% 
#  mutate(my = lubridate::as_date(paste(.$year, .$month, 01, sep = "-"))) %>% 
#  mutate(., cat = case_when(index >= 2 ~ "Extreme Wet",
#                            index < 2 & index >= 1.5 ~ "Severe Wet",
#                            index < 1.5 & index >= 1 ~ "Moderate Wet",
#                            index < 1 & index > -1 ~ "Normal",
#                            index <= -1 & index > -1.5 ~ "Moderate Drought",
#                            index <= -1.5 & index > -2 ~ "Severe Drought",
#                            index <= -2 ~ "Extreme Drought"),
#         type = "SPEI")
#gridmet.spei.df$cat<- factor(gridmet.spei.df$cat, levels = c("Extreme Wet","Severe Wet",
#                                                             "Moderate Wet","Normal",
#                                                             "Moderate Drought","Severe Drought",
#                                                             "Extreme Drought"))

#write.csv(gridmet.spei.df, "gridmet.spei_monthly.csv")

#3 categories, 1 as threshold:
#gridmet.spei.df.3cat<- gridmet.spei %>% 
#  map2(.x = ., .y = split(gridmet2, gridmet2$MU), ~{
#    .$fitted %>% 
#      data.frame() %>% 
#      rename(index = ".") %>% 
#      mutate(year = .y$year, month = .y$month)
#  }) %>% 
#  bind_rows(.id = "MU") %>% 
#  mutate(my = lubridate::as_date(paste(.$year, .$month, 01, sep = "-"))) %>%
#  filter(month >1,
#         month <=7) %>%
#  group_by(MU, year) %>% 
#  summarise(index = mean(index)) %>% 
#  mutate(., cat = case_when(index >= 1 ~ "Wet",
#                            index < 1 & index > -1 ~ "Normal",
#                            index <= -1 ~ "Drought"),
#         type = "SPEI")
#gridmet.spei.df.3cat$cat<- factor(gridmet.spei.df.3cat$cat, levels = c("Wet","Normal",
#                                                                       "Drought"))


#3 categories based on average values between Feb and Jul
#gridmet.spei.df.3cat_av02.07<- gridmet.spei %>% 
#  map2(.x = ., .y = split(gridmet2, gridmet2$MU), ~{
#    .$fitted %>% 
#      data.frame() %>% 
#      rename(index = ".") %>% 
#      mutate(year = .y$year, month = .y$month)
#  }) %>% 
#  bind_rows(.id = "MU") %>% 
#  mutate(my = lubridate::as_date(paste(.$year, .$month, 01, sep = "-"))) %>%
#  filter(month >1,
#         month <=7) %>%
#  group_by(MU, year) %>% 
#  summarise(index = mean(index)) %>% 
#  mutate(., cat = case_when(index >= 1 ~ "Wet",
#                            index < 1 & index > -1 ~ "Normal",
#                            index <= -1 ~ "Drought"),
#         type = "SPEI")
#gridmet.spei.df.3cat$cat<- factor(gridmet.spei.df.3cat$cat, levels = c("Wet","Normal",
#                                                             "Drought"))

#write.csv(gridmet.spei.df.3cat_av02.07, "spei_per_year_avg02-07_long.csv")

#3 categories, only extreme values (1.5 as threshold)
gridmet.spei.df.3cat_extreme <- gridmet.spei %>% 
  map2(.x = ., .y = split(gridmet2, gridmet2$MU), ~{
    .$fitted %>% 
      data.frame() %>% 
      rename(index = ".") %>% 
      mutate(year = .y$year, month = .y$month)
  }) %>% 
  bind_rows(.id = "MU") %>% 
  mutate(my = lubridate::as_date(paste(.$year, .$month, 01, sep = "-"))) %>% 
  mutate(., cat = case_when(index >= 1.5 ~ "Wet",
                            index < 1.5 & index > -1.5 ~ "Normal",
                            index <= -1.5 ~ "Drought"),
         type = "SPEI")
gridmet.spei.df.3cat_extreme$cat <- factor(gridmet.spei.df.3cat_extreme$cat, 
                                           levels = c("Wet","Normal","Drought"))



# Plot SPEI/SPI time series
ggplot(gridmet.spei.df.3cat_extreme %>% 
         filter(month == 5), aes(year, index, color = MU)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = c(-1.5,1.5), linetype = "dotted", color = "red") +
  geom_line(size = 1) +
  scale_color_manual("", values = pal) +
  scale_x_continuous(breaks = c(1995, 2000, 2005, 2010, 2015, 2020, 2025)) +
  labs(x = "Year", y = "SP(E)I", title = "18 Month (Dec - May) Temporal Period",
       subtitle = "Data source: gridMET") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14, face = "bold"),
        strip.background = element_blank(),
        legend.text = element_text(size = 12)) 

ggsave(filename = "figures/StaPre_EvaInd_5_Regions.jpg", dpi = 300)

ggplot(gridmet.spei.df.3cat_extreme %>% 
         filter(month == 5), aes(year, index, color = type)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(size = 1) +
  # scale_color_manual("", values = pal) +
  scale_color_brewer("", palette = "Dark2") +
  labs(x = "Year", y = "SP(E)I", title = "18 Month (Dec - May) Temporal Period",
       subtitle = "Data source: gridMET") +
  theme_bw() +
#  xlim(1996,2025) +
  theme(#panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14, face = "bold"),
        strip.background = element_blank(),
        legend.text = element_text(size = 12)) +
  facet_wrap(~MU)


ggplot(gridmet.spei.df.3cat_extreme %>% 
         filter(month == 5, year > 1995), aes(year, index)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(size = 1) +
  geom_point(aes(fill = cat), size = 3, shape = 21) +
  scale_fill_brewer("", palette = 'RdBu', direction = -1) +
  scale_x_continuous(breaks = c(1995, 2000, 2005, 2010, 2015, 2020, 2025)) +
  labs(x = "Year", y = "SP(E)I", title = "18 Month (Dec - May) Temporal Period",
       subtitle = "Data source: gridMET") +
  theme_bw() +
#  xlim(1996, 2024) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14, face = "bold"),
        strip.background = element_blank(),
        legend.text = element_text(size = 12)) +
  facet_grid(MU ~ type)


# SPEI per year  ----------------------------------------------------------
#May values
#threshold =1
#to visualize
#gridmet.spei.df.yr <- gridmet.spei.df %>% 
#  filter(month == 5) %>% 
#  select(MU,
#         index,
#         year) %>% 
#  pivot_wider(values_from = index,
#              names_from = MU)

#write.csv(gridmet.spei.df.yr, "gridmet.spei.yr.csv")

#table in long format for modeling
#gridmet.spei.df.yr <-  gridmet.spei.df.3cat %>% 
#  filter (month == 5) 

#using average values Feb-Jul
#gridmet.spei.df.yr_av02.07 <- gridmet.spei.df.3cat_av02.07 %>% 
#  select(MU,
#         cat,
#         year) %>% 
#  pivot_wider(values_from = cat,
#              names_from = MU)


#write.csv(gridmet.spei.df.yr_av02.07, "spei_per_year_avg02-07.csv")

#table in long format for modeling
#gridmet.spei.df.yr_long <- gridmet.spei.df.3cat %>% 
#  group_by(MU, year) %>% 
#  summarise(index = mean(index)) %>% 
#  select(MU,
#         index,
#         year)


#gridmet.spei.df.yr.table <- gridmet.spei.df.yr %>% 
#  select(MU,
#         index,
#         year,
#         cat) 


#write.csv(gridmet.spei.df.yr.table, "SPEI_per_yr_long.csv")

#threshold: 1.5
gridmet.spei.df.yr_ex <-  gridmet.spei.df.3cat_extreme %>% 
  filter (month == 5) 

gridmet.spei.df.yr.table_ex <- gridmet.spei.df.yr_ex %>% 
  dplyr::select(MU,
         year,
         cat) 



## add invasion data -------------------------------------------------------
inv <- read.csv("data/tsi_by_region_for_ss_96_23_50perc.csv")
inv$invaded <- ifelse(inv$TSI < 0, 0, 1) 
inv.data <- inv %>% 
  mutate(inv.cat = case_when(invaded == 0 ~ "No",
                             invaded == 1 ~ "Yes"))

head(inv.data)

# This data is up to 2023, but the last year of invasion is 2018, so between 2018 and 2023 all regions has been invaded
inv.data |>
  filter(Year %in% c(2018:2023)) |>
  group_by(Region) |>
  summarise(inv.cat = unique(inv.cat))

# Can we assume that in 2024-2025 the regions has continued being invaded?


gridmet.spei.df.yr.table_ex_inv <- gridmet.spei.df.yr.table_ex %>%  
  left_join(inv.data, join_by(year==Year, MU == Region)) %>%
  rename(spei=cat,
         Pop = MU) %>% 
  dplyr::select(Pop,
                spei,
                inv.cat,
                year)

head(gridmet.spei.df.yr.table_ex_inv)

write.csv(gridmet.spei.df.yr.table_ex_inv, "data/SPEI_per_yr_long_extreme_inv.csv")

#for visualization
gridmet.spei.df.yr.table_ex_wide <- gridmet.spei.df.yr.table_ex %>% 
  pivot_wider(values_from = cat,
              names_from = MU)


inv.data_wide <- inv.data %>% 
  dplyr::select(Year,
                Region,
                inv.cat) %>% 
  pivot_wider(values_from = inv.cat,
              names_from = Region)

gridmet.spei.df.yr.table_ex_wide_inv <- gridmet.spei.df.yr.table_ex_wide %>% left_join(inv.data_wide, join_by(year==Year))



#contrast with original table
#Table2_original <- read.csv("Table2_original.csv")


#gridmet.spei.df.yr.table2_2020 <- gridmet.spei.df.yr.table2 %>% 
#  filter(year < 2021 & year >1995)

#which(gridmet.spei.df.yr.table2_2020!= Table2_original, arr.ind =TRUE)


# check current data ------------------------------------------------------

gridmet.spei.df.yr.2020 <- gridmet.spei.df.yr.table %>% 
  filter (year > 2020,
          cat == "Drought")

#-1.165 and -1.131

#plot
ggplot(gridmet.spei.df %>% 
         filter(month == 5), aes(year, index, color = MU)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(size = 1) +
  scale_color_manual("", values = pal) +
  scale_x_continuous(breaks = c(1995, 2000, 2005, 2010, 2015, 2020, 2025)) +
  labs(x = "Year", y = "SP(E)I", title = "18 Month (Dec - May) Temporal Period",
       subtitle = "Data source: gridMET") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14, face = "bold"),
        strip.background = element_blank(),
        legend.text = element_text(size = 12)) 


gridmet.spei.df.yr.Mod_d <- gridmet.spei.df %>% 
  filter (month == 5,
          index <= -1 & index > -1.5)

#13 sites-years with moderate drought

gridmet.spei.df.yr.Mod_w <- gridmet.spei.df %>% 
  filter (month == 5,
          index >= 1 & index < 1.5)

#17 sites-years with moderate wet




