
#Download MACA data (from all available GCMs) using 'climateR' package
MACA_summary = function(AOI, params, scenarios, startDate, endDate) {
  # AOI = an sf or sp object on which to extract climate data
  # params = a vector of the model params to download from MACA
  # scenarios = vector of RCP climate scenarios to evaluate
  # startDate = a character string of the start date in ISO format
  # endDate = a character string of the end date in ISO format
  # timeRes = a character of 'daily' or 'monthly' to download data
  
  # tictoc::tic()
  ensembles = getMACA(AOI = AOI,
                      param = params, 
                      model = model_meta$maca$model,
                      scenario = scenarios,
                      startDate = startDate, 
                      endDate = endDate,
                      timeRes = "daily")
  # tictoc::toc()
  print(paste("Finished downloading data for", AOI$MU_Name))
  
  
  # Calculate sum total rainfall (mm) for precip and avg min/max temps
  if (lubridate::year(as.Date(startDate, format = "%Y-%m-%d")) <= 2005) {
    scenarios<- "historical"
  }
  
  clim<- vector('list', length = length(params) * length(scenarios))
  oo<- 1  #create counter to store results in list
  for (i in 1:length(params)) { #print(paste("param",i))
    for (j in 1:length(scenarios)) { #print(paste("scenario",j))
      
      # index combos of params and scenarios
      if (lubridate::year(as.Date(startDate, format = "%Y-%m-%d")) <= 2005) {
        name1<- paste(params[i], "historical", sep = "_")  #need to index historical data differently
      } else {
        name1<- paste(params[i], scenarios[j], sep = "_")
      }
      
      ind<- grep(name1, names(ensembles))
      
      # perform function (sum or mean) on each set of predictions over time period
      if (params[i] == 'prcp') {
        clim[[oo]]<- purrr::map(ensembles[ind], ~daily2monthly(., fun = sum)) %>% 
          raster::brick()
        names(clim)[oo]<- name1
        
      } else {
        clim[[oo]]<- purrr::map(ensembles[ind], ~daily2monthly(., fun = mean)) %>% 
          raster::brick()
        names(clim)[oo]<- name1
        
      }  #close if-else conditional
      
      oo<- oo + 1
    }  # close j
  }  #close i
  
  
  # create vector of model names
  n.yrmonths<- nlayers(clim[[1]])/20
  gcm<- sort(model_meta$maca$model[1:20]) %>% 
    rep(., each = n.yrmonths)
  
  # start.month<- lubridate::month(as.Date(startDate, format = "%Y-%m-%d"))
  # end.month<- lubridate::month(as.Date(endDate, format = "%Y-%m-%d"))
  date.seq<- seq.Date(as.Date(startDate, format = "%Y-%m-%d"),
                      as.Date(endDate, format = "%Y-%m-%d"),
                      by = "month")
  months1<- rep(lubridate::month(date.seq), 20)
  year1<- rep(lubridate::year(date.seq), 20) 
  
  
  if (lubridate::year(as.Date(startDate, format = "%Y-%m-%d")) <= 2005) {
    # summarize params by model for historical
    tmp.prcp<- t(exact_extract(clim$prcp_historical, AOI, "mean"))
    tmp.tmin<- t(exact_extract(clim$tmin_historical, AOI, "mean"))
    tmp.tmax<- t(exact_extract(clim$tmax_historical, AOI, "mean"))
    mod.res<- data.frame(model = gcm, month = months1, year = year1, prcp = tmp.prcp,
                         tmin = tmp.tmin, tmax = tmp.tmax, scenario = "historical",
                         row.names = NULL)
    
  } else {
    # summarize params by model for RCP 4.5
    tmp.prcp<- t(exact_extract(clim$prcp_rcp45, AOI, "mean"))
    tmp.tmin<- t(exact_extract(clim$tmin_rcp45, AOI, "mean"))
    tmp.tmax<- t(exact_extract(clim$tmax_rcp45, AOI, "mean"))
    mod.res45<- data.frame(model = gcm, month = months1, year = year1, prcp = tmp.prcp,
                           tmin = tmp.tmin, tmax = tmp.tmax, scenario = "rcp45", row.names = NULL)
    
    # summarize params by model for RCP 8.5
    tmp.prcp<- t(exact_extract(clim$prcp_rcp85, AOI, "mean"))
    tmp.tmin<- t(exact_extract(clim$tmin_rcp85, AOI, "mean"))
    tmp.tmax<- t(exact_extract(clim$tmax_rcp85, AOI, "mean"))
    mod.res85<- data.frame(model = gcm, month = months1, year = year1, prcp = tmp.prcp,
                           tmin = tmp.tmin, tmax = tmp.tmax, scenario = "rcp85", row.names = NULL)
    
    
    # Combine all summary results in single DF
    mod.res<- rbind(mod.res45, mod.res85)
  }
  
  return(mod.res)
}

#----------------------------

# Internal function to change the time scale of the environ covars
daily2monthly = function(stack, fun) {
  # stack = a RasterStack on which to convert var from daily to monthly scale (per year)
  # fun = a function specifying how to aggregate the daily data to monthly scale
  
  #get the date from the names of the layers and extract the year-month
  ind <- format(as.Date(names(stack), format = "X%Y.%m.%d"), format = "%Y-%m")
  # ind <- as.numeric(ind)
  
  #apply function on layers
  tmp<- stackApply(stack, ind, fun = fun)
  names(tmp)<- gsub("index_", "X", names(tmp))
  
  return(tmp)
}

#----------------------------

#Download gridMET data using the 'cimateR' package
gridMET_summary = function(AOI, params, startDate, endDate) {
  # AOI = an sf or sp object on which to extract climate data
  # params = a vector of the model params to download from MACA
  # startDate = a character string of the start date in ISO format
  # endDate = a character string of the end date in ISO format
  
  # tictoc::tic()
  gridmet = getGridMET(AOI = AOI,
                       param = params, 
                       startDate = startDate, 
                       endDate = endDate)
  # tictoc::toc()
  print(paste("Finished downloading data for", AOI$MU_Name))
  
  
  # Calculate sum total rainfall (mm) for precip and avg min/max temps
  clim<- vector('list', length = length(params))
  for (i in 1:length(params)) { #print(paste("param",i))
    ind<- grep(params[i], names(gridmet))
    
    # perform function (sum or mean) on each set of predictions over time period
    if (params[i] == 'prcp') {
      clim[[i]]<- purrr::map(gridmet[ind], ~daily2monthly(., fun = sum)) %>% 
        raster::brick()
      names(clim)[i]<- params[i]
      
    } else if (params[i] == 'palmer') {
      clim[[i]]<- purrr::map(gridmet[ind], ~pentad2monthly(., fun = mean)) %>% 
        raster::brick()
      names(clim)[i]<- params[i]
      
    } else {
      clim[[i]]<- purrr::map(gridmet[ind], ~daily2monthly(., fun = mean)) %>% 
        raster::brick()
      names(clim)[i]<- params[i]
      
    }  #close if-else conditional
  }  #close i
  
  
  # create vector of model names
  n.yrmonths<- nlayers(clim[[1]])
  
  date.seq<- seq.Date(as.Date(startDate, format = "%Y-%m-%d"),
                      as.Date(endDate, format = "%Y-%m-%d"),
                      by = "month")
  months1<- lubridate::month(date.seq)
  year1<- lubridate::year(date.seq)
  
  
  # summarize params for gridMET
  if (!("palmer" %in% params)) {
    tmp.prcp<- t(exact_extract(clim$prcp, AOI, "mean"))
    tmp.tmin<- t(exact_extract(clim$tmin, AOI, "mean"))
    tmp.tmax<- t(exact_extract(clim$tmax, AOI, "mean"))
    clim.res<- data.frame(month = months1, year = year1, prcp = tmp.prcp,
                          tmin = tmp.tmin, tmax = tmp.tmax, row.names = NULL)
    
  } else if ("palmer" %in% params) {
    tmp.palmer<- t(exact_extract(clim$palmer, AOI, "mean"))
    clim.res<- data.frame(month = months1, year = year1, pdsi = tmp.palmer,
                          row.names = NULL)
  }
  
  
  return(clim.res)
}

#----------------------------

# Internal function for dealing with 10-day PDSI data
pentad2monthly = function(stack, fun) {
  # stack = a RasterStack on which to convert var from 10-day to monthly scale (per year)
  # fun = a function specifying how to aggregate the 10-day data to monthly scale
  
  #get the date from the names of the layers and extract the year-month
  ind<- names(stack)
  
  for (i in 1:nlayers(stack)) {
    month1<- dplyr::case_when(i %in% 1:6 ~ "01",
                              i %in% 7:12 ~ "02",
                              i %in% 13:18 ~ "03",
                              i %in% 19:24 ~ "04",
                              i %in% 25:30 ~ "05",
                              i %in% 31:36 ~ "06",
                              i %in% 37:42 ~ "07",
                              i %in% 43:48 ~ "08",
                              i %in% 49:54 ~ "09",
                              i %in% 55:60 ~ "10",
                              i %in% 61:66 ~ "11",
                              i %in% 67:73 ~ "12")
    ind<- gsub(paste0("pentad_",i,"$"), month1, ind)
  }
  
  
  #apply function on layers
  tmp<- stackApply(stack, ind, fun = fun)
  names(tmp)<- gsub("index_", "", names(tmp))
  
  return(tmp)
}


#----------------------------

# purrr::flatten, but with ability to combine names from nested list
flatten2 <- function(x, keep_names = TRUE, sep = "_") {
  
  # Handling unnamed list & keep_names = TRUE
  if (keep_names & all(is.null(names(x)))) {
    abort("You can't keep names of an unnamed list")
  }
  
  res <- purrr::flatten(x)
  
  if (keep_names)  {
    rep1<- dplyr::n_distinct(names(res))
    nm <- paste(rep(names(x), each = rep1), names(res), sep = sep)
    return(purrr::set_names(res, nm))
  }
  res
}

#---------------------------

# Calculates B matrix using survival and repro rates
getB <- function(x, spei){
  # x = a data frame with columns for MU, age_cat, spei_cat, and demographic params
  # spei = a character string denoting the SPEI category from which to draw demographic rates
  
  
  B0 <- matrix(0, Ns, Ns)  						             # empty blocks
  
  #filter by spei
  tmp<- x %>% 
    filter(spei_cat == spei)
  
  #filter by age class
  tmp.a<- tmp %>% 
    filter(age_cat == "adult")
  tmp.j<- tmp %>% 
    filter(age_cat == "juv")
  
  #main code convert 0 to 1e-5 to get reliable values from beta (won't take 0)
  
  #define params
  Sj <- getBetaDist(meanx = tmp.j$S, stdevx = tmp.j$S.SE)
  Sa <- getBetaDist(meanx = tmp.a$S, stdevx = tmp.a$S.SE)
  NS <- getBetaDist(meanx = tmp.a$NS, stdevx = tmp.a$NS_SE)
  Jbreed <- getBetaDist(meanx = tmp.j$Breed, stdevx = tmp.j$Breed.SE)
  Abreed <- getBetaDist(meanx = tmp.a$Breed, stdevx = tmp.a$Breed.SE)
  Jattempt <- getlogNormDist(meanx = tmp.j$Attempt, stdevx = tmp.j$Attempt.SE)
  Aattempt <- getlogNormDist(meanx = tmp.a$Attempt, stdevx = tmp.a$Attempt.SE)
  Y <- getlogNormDist(meanx = tmp.a$Yng_cnt, stdevx = tmp.a$Yng_cnt.SE)
  
  #adjust very small demographic rates to be zero
  Sj[Sj < 1e-5] <- 0
  Sa[Sa < 1e-5] <- 0
  NS[NS < 1e-5] <- 0
  Jbreed[Jbreed < 1e-5] <- 0
  Abreed[Abreed < 1e-5] <- 0
  
  #adjust breeding probability for carrying capacity
  # ns.t.mat <- matrix(N.t.i[[k]][,i-1], nrow=Np, ncol=Ns, byrow = T)
  # ns.t.total <- rowSums(ns.t.mat)
  # extra <- ns.t.total - K
  # Jbreed[extra > 0] <- 0 #for sites > K, no breeding
  # Abreed[extra > 0] <- 0 #for sites > K, no breeding
  
  #fecundity
  ms <- Jbreed * NS * Y/2 * Jattempt
  ma <- Abreed * NS * Y/2 * Aattempt
  
  #fertility
  Fs <- Sj * ms
  Fa <- Sj * ma
  
  EAST <- matrix(c(Fs[1], Fa[1], Sa[1], Sa[1]), Ns, Ns, byrow=T)
  EVER <- matrix(c(Fs[2], Fa[2], Sa[2], Sa[2]), Ns, Ns, byrow=T)
  KRV <-  matrix(c(Fs[3], Fa[3], Sa[3], Sa[3]), Ns, Ns, byrow=T)
  OKEE <- matrix(c(Fs[4], Fa[4], Sa[4], Sa[4]), Ns, Ns, byrow=T)
  PP <- matrix(c(Fs[5], Fa[5], Sa[5], Sa[5]), Ns, Ns, byrow=T)
  SJM <- matrix(c(Fs[6], Fa[6], Sa[6], Sa[6]), Ns, Ns, byrow=T)
  
  #combine: number of blocks = number of sites
  B <- rbind(cbind(EAST,B0,B0,B0,B0,B0),
             cbind(B0,EVER,B0,B0,B0,B0),
             cbind(B0,B0,KRV,B0,B0,B0),
             cbind(B0,B0,B0,OKEE,B0,B0),
             cbind(B0,B0,B0,B0,PP,B0),
             cbind(B0,B0,B0,B0,B0,SJM))
  return(B)
}


# Calculates B matrix using survival and repro rates

getB_past <- function(x, spei){
  # x = a data frame with columns for MU, age_cat, spei_cat, and demographic params
  # spei = a vector with character strings denoting the SPEI category at each MU from which to draw demographic rates
  
  
  B0 <- matrix(0, Ns, Ns)  						             # empty blocks
  

  #filter by age class
  tmp.a<- x %>% 
    filter(age_cat == "adult")
  tmp.j<- x %>% 
    filter(age_cat == "juv")
  
  #main code convert 0 to 1e-5 to get reliable values from beta (won't take 0)
  
  #define params
  Sj <- getBetaDist(meanx = tmp.j$S, stdevx = tmp.j$S.SE)
  Sa <- getBetaDist(meanx = tmp.a$S, stdevx = tmp.a$S.SE)
  NS <- getBetaDist(meanx = tmp.a$NS, stdevx = tmp.a$NS_SE)
  Jbreed <- getBetaDist(meanx = tmp.j$Breed, stdevx = tmp.j$Breed.SE)
  Abreed <- getBetaDist(meanx = tmp.a$Breed, stdevx = tmp.a$Breed.SE)
  Jattempt <- getlogNormDist(meanx = tmp.j$Attempt, stdevx = tmp.j$Attempt.SE)
  Aattempt <- getlogNormDist(meanx = tmp.a$Attempt, stdevx = tmp.a$Attempt.SE)
  Y <- getlogNormDist(meanx = tmp.a$Yng_cnt, stdevx = tmp.a$Yng_cnt.SE)
  
  #adjust very small demographic rates to be zero
  Sj[Sj < 1e-5] <- 0
  Sa[Sa < 1e-5] <- 0
  NS[NS < 1e-5] <- 0
  Jbreed[Jbreed < 1e-5] <- 0
  Abreed[Abreed < 1e-5] <- 0
  
  #adjust breeding probability for carrying capacity
  # ns.t.mat <- matrix(N.t.i[[k]][,i-1], nrow=Np, ncol=Ns, byrow = T)
  # ns.t.total <- rowSums(ns.t.mat)
  # extra <- ns.t.total - K
  # Jbreed[extra > 0] <- 0 #for sites > K, no breeding
  # Abreed[extra > 0] <- 0 #for sites > K, no breeding
  
  #fecundity
  ms <- Jbreed * NS * Y/2 * Jattempt
  ma <- Abreed * NS * Y/2 * Aattempt
  
  #fertility
  Fs <- Sj * ms
  Fa <- Sj * ma
  
  if(spei[1] == "Drought") {
    EAST <- matrix(c(Fs[1], Fa[1], Sa[1], Sa[1]), Ns, Ns, byrow=T)
  }
  if(spei[1] == "Normal") {
    EAST <- matrix(c(Fs[2], Fa[2], Sa[2], Sa[2]), Ns, Ns, byrow=T)
  }
  if(spei[1] == "Wet") {
    EAST <- matrix(c(Fs[3], Fa[3], Sa[3], Sa[3]), Ns, Ns, byrow=T)
  }
  
  
  if(spei[2] == "Drought") {
    EVER <- matrix(c(Fs[4], Fa[4], Sa[4], Sa[4]), Ns, Ns, byrow=T)
  }
  if(spei[2] == "Normal") {
    EVER <- matrix(c(Fs[5], Fa[5], Sa[5], Sa[5]), Ns, Ns, byrow=T)
  }
  if(spei[2] == "Wet") {
    EVER <- matrix(c(Fs[6], Fa[6], Sa[6], Sa[6]), Ns, Ns, byrow=T)
  }
  
  if(spei[3] == "Drought") {
    KRV <- matrix(c(Fs[7], Fa[7], Sa[7], Sa[7]), Ns, Ns, byrow=T)
  }
  if(spei[3] == "Normal") {
    KRV <- matrix(c(Fs[8], Fa[8], Sa[8], Sa[8]), Ns, Ns, byrow=T)
  }
  if(spei[3] == "Wet") {
    KRV <- matrix(c(Fs[9], Fa[9], Sa[9], Sa[9]), Ns, Ns, byrow=T)
  }
  
  if(spei[4] == "Drought") {
    OKEE <- matrix(c(Fs[10], Fa[10], Sa[10], Sa[10]), Ns, Ns, byrow=T)
  }
  if(spei[4] == "Normal") {
    OKEE <- matrix(c(Fs[11], Fa[11], Sa[11], Sa[11]), Ns, Ns, byrow=T)
  }
  if(spei[4] == "Wet") {
    OKEE <- matrix(c(Fs[12], Fa[12], Sa[12], Sa[12]), Ns, Ns, byrow=T)
  }
  
  if(spei[5] == "Drought") {
    PP <- matrix(c(Fs[13], Fa[13], Sa[13], Sa[13]), Ns, Ns, byrow=T)
  }
  if(spei[5] == "Normal") {
    PP <- matrix(c(Fs[14], Fa[14], Sa[14], Sa[14]), Ns, Ns, byrow=T)
  }
  if(spei[5] == "Wet") {
    PP <- matrix(c(Fs[15], Fa[15], Sa[15], Sa[15]), Ns, Ns, byrow=T)
  }
  
  if(spei[6] == "Drought") {
    SJM <- matrix(c(Fs[16], Fa[16], Sa[16], Sa[16]), Ns, Ns, byrow=T)
  }
  if(spei[6] == "Normal") {
    SJM <- matrix(c(Fs[17], Fa[17], Sa[17], Sa[17]), Ns, Ns, byrow=T)
  }
  if(spei[6] == "Wet") {
    SJM <- matrix(c(Fs[18], Fa[18], Sa[18], Sa[18]), Ns, Ns, byrow=T)
  }
  
  
  
  #combine: number of blocks = number of sites
  B <- rbind(cbind(EAST,B0,B0,B0,B0,B0),
             cbind(B0,EVER,B0,B0,B0,B0),
             cbind(B0,B0,KRV,B0,B0,B0),
             cbind(B0,B0,B0,OKEE,B0,B0),
             cbind(B0,B0,B0,B0,PP,B0),
             cbind(B0,B0,B0,B0,B0,SJM))
  return(B)
}



getB_floater_past <- function(x, spei){
  # x = a data frame with columns for MU, age_cat, spei_cat, and demographic params
  # spei = a character string denoting the SPEI category from which to draw demographic rates
  
  
  B0 <- matrix(0, Ns, Ns)  						             # empty blocks
  
  #filter by age class
  tmp.a<- x %>% 
    filter(age_cat == "adult")
  tmp.j<- x %>% 
    filter(age_cat == "juv")
  
  #main code convert 0 to 1e-5 to get reliable values from beta (won't take 0)
  
  #define params
  Sj <- getBetaDist(meanx = tmp.j$S, stdevx = tmp.j$S.SE)
  Sa <- getBetaDist(meanx = tmp.a$S, stdevx = tmp.a$S.SE)
  NS <- getBetaDist(meanx = tmp.a$NS, stdevx = tmp.a$NS_SE)
  Jbreed <- getBetaDist(meanx = tmp.j$Breed, stdevx = tmp.j$Breed.SE)
  Abreed <- getBetaDist(meanx = tmp.a$Breed, stdevx = tmp.a$Breed.SE)
  Jattempt <- getlogNormDist(meanx = tmp.j$Attempt, stdevx = tmp.j$Attempt.SE)
  Aattempt <- getlogNormDist(meanx = tmp.a$Attempt, stdevx = tmp.a$Attempt.SE)
  Y <- getlogNormDist(meanx = tmp.a$Yng_cnt, stdevx = tmp.a$Yng_cnt.SE)
  
  #adjust very small demographic rates to be zero
  Sj[Sj < 1e-5] <- 0
  Sa[Sa < 1e-5] <- 0
  NS[NS < 1e-5] <- 0
  Jbreed[Jbreed < 1e-5] <- 0
  Abreed[Abreed < 1e-5] <- 0
  
  #adjust breeding probability for carrying capacity
  # ns.t.mat <- matrix(N.t.i[[k]][,i-1], nrow=Np, ncol=Ns, byrow = T)
  # ns.t.total <- rowSums(ns.t.mat)
  # extra <- ns.t.total - K
  # Jbreed[extra > 0] <- 0 #for sites > K, no breeding
  # Abreed[extra > 0] <- 0 #for sites > K, no breeding
  
  #fecundity: no breeding for floaters 
  ms <- 0 * NS * Y/2 * Jattempt 
  ma <- 0 * NS * Y/2 * Aattempt 
  
  #fertility
  Fs <- Sj * ms
  Fa <- Sj * ma
  
  if(spei[1] == "Drought") {
    EAST <- matrix(c(Fs[1], Fa[1], Sa[1], Sa[1]), Ns, Ns, byrow=T)
  }
  if(spei[1] == "Normal") {
    EAST <- matrix(c(Fs[2], Fa[2], Sa[2], Sa[2]), Ns, Ns, byrow=T)
  }
  if(spei[1] == "Wet") {
    EAST <- matrix(c(Fs[3], Fa[3], Sa[3], Sa[3]), Ns, Ns, byrow=T)
  }
  
  
  if(spei[2] == "Drought") {
    EVER <- matrix(c(Fs[4], Fa[4], Sa[4], Sa[4]), Ns, Ns, byrow=T)
  }
  if(spei[2] == "Normal") {
    EVER <- matrix(c(Fs[5], Fa[5], Sa[5], Sa[5]), Ns, Ns, byrow=T)
  }
  if(spei[2] == "Wet") {
    EVER <- matrix(c(Fs[6], Fa[6], Sa[6], Sa[6]), Ns, Ns, byrow=T)
  }
  
  if(spei[3] == "Drought") {
    KRV <- matrix(c(Fs[7], Fa[7], Sa[7], Sa[7]), Ns, Ns, byrow=T)
  }
  if(spei[3] == "Normal") {
    KRV <- matrix(c(Fs[8], Fa[8], Sa[8], Sa[8]), Ns, Ns, byrow=T)
  }
  if(spei[3] == "Wet") {
    KRV <- matrix(c(Fs[9], Fa[9], Sa[9], Sa[9]), Ns, Ns, byrow=T)
  }
  
  if(spei[4] == "Drought") {
    OKEE <- matrix(c(Fs[10], Fa[10], Sa[10], Sa[10]), Ns, Ns, byrow=T)
  }
  if(spei[4] == "Normal") {
    OKEE <- matrix(c(Fs[11], Fa[11], Sa[11], Sa[11]), Ns, Ns, byrow=T)
  }
  if(spei[4] == "Wet") {
    OKEE <- matrix(c(Fs[12], Fa[12], Sa[12], Sa[12]), Ns, Ns, byrow=T)
  }
  
  if(spei[5] == "Drought") {
    PP <- matrix(c(Fs[13], Fa[13], Sa[13], Sa[13]), Ns, Ns, byrow=T)
  }
  if(spei[5] == "Normal") {
    PP <- matrix(c(Fs[14], Fa[14], Sa[14], Sa[14]), Ns, Ns, byrow=T)
  }
  if(spei[5] == "Wet") {
    PP <- matrix(c(Fs[15], Fa[15], Sa[15], Sa[15]), Ns, Ns, byrow=T)
  }
  
  if(spei[6] == "Drought") {
    SJM <- matrix(c(Fs[16], Fa[16], Sa[16], Sa[16]), Ns, Ns, byrow=T)
  }
  if(spei[6] == "Normal") {
    SJM <- matrix(c(Fs[17], Fa[17], Sa[17], Sa[17]), Ns, Ns, byrow=T)
  }
  if(spei[6] == "Wet") {
    SJM <- matrix(c(Fs[18], Fa[18], Sa[18], Sa[18]), Ns, Ns, byrow=T)
  }
  
  
  #combine: number of blocks = number of sites
  B <- rbind(cbind(EAST,B0,B0,B0,B0,B0),
             cbind(B0,EVER,B0,B0,B0,B0),
             cbind(B0,B0,KRV,B0,B0,B0),
             cbind(B0,B0,B0,OKEE,B0,B0),
             cbind(B0,B0,B0,B0,PP,B0),
             cbind(B0,B0,B0,B0,B0,SJM))
  return(B)
}


getB_floater <- function(x, spei){
  # x = a data frame with columns for MU, age_cat, spei_cat, and demographic params
  # spei = a character string denoting the SPEI category from which to draw demographic rates
  
  
  B0 <- matrix(0, Ns, Ns)  						             # empty blocks
  
  #filter by spei
  tmp<- x %>% 
    filter(spei_cat == spei)
  
  #filter by age class
  tmp.a<- tmp %>% 
    filter(age_cat == "adult")
  tmp.j<- tmp %>% 
    filter(age_cat == "juv")
  
  #main code convert 0 to 1e-5 to get reliable values from beta (won't take 0)
  
  #define params
  Sj <- getBetaDist(meanx = tmp.j$S, stdevx = tmp.j$S.SE)
  Sa <- getBetaDist(meanx = tmp.a$S, stdevx = tmp.a$S.SE)
  NS <- getBetaDist(meanx = tmp.a$NS, stdevx = tmp.a$NS_SE)
  Jbreed <- getBetaDist(meanx = tmp.j$Breed, stdevx = tmp.j$Breed.SE)
  Abreed <- getBetaDist(meanx = tmp.a$Breed, stdevx = tmp.a$Breed.SE)
  Jattempt <- getlogNormDist(meanx = tmp.j$Attempt, stdevx = tmp.j$Attempt.SE)
  Aattempt <- getlogNormDist(meanx = tmp.a$Attempt, stdevx = tmp.a$Attempt.SE)
  Y <- getlogNormDist(meanx = tmp.a$Yng_cnt, stdevx = tmp.a$Yng_cnt.SE)
  
  #adjust very small demographic rates to be zero
  Sj[Sj < 1e-5] <- 0
  Sa[Sa < 1e-5] <- 0
  NS[NS < 1e-5] <- 0
  Jbreed[Jbreed < 1e-5] <- 0
  Abreed[Abreed < 1e-5] <- 0
  
  #adjust breeding probability for carrying capacity
  # ns.t.mat <- matrix(N.t.i[[k]][,i-1], nrow=Np, ncol=Ns, byrow = T)
  # ns.t.total <- rowSums(ns.t.mat)
  # extra <- ns.t.total - K
  # Jbreed[extra > 0] <- 0 #for sites > K, no breeding
  # Abreed[extra > 0] <- 0 #for sites > K, no breeding
  
  #fecundity: no breeding for floaters 
  ms <- 0 * NS * Y/2 * Jattempt 
  ma <- 0 * NS * Y/2 * Aattempt 
  
  #fertility
  Fs <- Sj * ms
  Fa <- Sj * ma
  
  EAST <- matrix(c(Fs[1], Fa[1], Sa[1], Sa[1]), Ns, Ns, byrow=T)
  EVER <- matrix(c(Fs[2], Fa[2], Sa[2], Sa[2]), Ns, Ns, byrow=T)
  KRV <-  matrix(c(Fs[3], Fa[3], Sa[3], Sa[3]), Ns, Ns, byrow=T)
  OKEE <- matrix(c(Fs[4], Fa[4], Sa[4], Sa[4]), Ns, Ns, byrow=T)
  PP <- matrix(c(Fs[5], Fa[5], Sa[5], Sa[5]), Ns, Ns, byrow=T)
  SJM <- matrix(c(Fs[6], Fa[6], Sa[6], Sa[6]), Ns, Ns, byrow=T)
  
  #combine: number of blocks = number of sites
  B <- rbind(cbind(EAST,B0,B0,B0,B0,B0),
             cbind(B0,EVER,B0,B0,B0,B0),
             cbind(B0,B0,KRV,B0,B0,B0),
             cbind(B0,B0,B0,OKEE,B0,B0),
             cbind(B0,B0,B0,B0,PP,B0),
             cbind(B0,B0,B0,B0,B0,SJM))
  return(B)
}

# Calculates B matrix using survival and repro rates
#only difference is the data passed for hindcasting
getB_hindcast <- function(x, MUspei){
  # x = a data frame with columns for MU, age_cat, spei_cat, and demographic params
  # MUspei = dataframe with SPEI category for each MU from which to draw demographic rates
  
  
  B0 <- matrix(0, Ns, Ns)  						             # empty blocks
  
  #filter by MUspei
  tmp<- merge(MUspei, x, by = c("MU", "spei_cat"))
  
  #filter by age class
  tmp.a<- tmp %>% 
    filter(age_cat == "adult")
  tmp.j<- tmp %>% 
    filter(age_cat == "juv")
  
  #main code convert 0 to 1e-5 to get reliable values from beta (won't take 0)
  
  #define params
  Sj <- getBetaDist(meanx = tmp.j$S, stdevx = tmp.j$S - tmp.j$S_LCL)
  Sa <- getBetaDist(meanx = tmp.a$S, stdevx = tmp.a$S - tmp.a$S_LCL)
  NS <- getBetaDist(meanx = tmp.a$NS, stdevx = tmp.a$NS_UCL - tmp.a$NS)
  Jbreed <- getBetaDist(meanx = tmp.j$Breed, stdevx = tmp.j$Breed.SE)
  Abreed <- getBetaDist(meanx = tmp.a$Breed, stdevx = tmp.a$Breed.SE)
  Jattempt <- getlogNormDist(meanx = tmp.j$Attempt, stdevx = tmp.j$Attempt.SE)
  Aattempt <- getlogNormDist(meanx = tmp.a$Attempt, stdevx = tmp.a$Attempt.SE)
  Y <- getlogNormDist(meanx = tmp.a$Yng_cnt, stdevx = tmp.a$Yng_cnt.SE)
  
  #adjust very small demographic rates to be zero
  Sj[Sj < 1e-5] <- 0
  Sa[Sa < 1e-5] <- 0
  NS[NS < 1e-5] <- 0
  Jbreed[Jbreed < 1e-5] <- 0
  Abreed[Abreed < 1e-5] <- 0
  
  #fecundity
  ms <- Jbreed * NS * Y/2 * Jattempt
  ma <- Abreed * NS * Y/2 * Aattempt
  
  #fertility
  Fs <- Sj * ms
  Fa <- Sj * ma
  
  EAST <- matrix(c(Fs[1], Fa[1], Sa[1], Sa[1]), Ns, Ns, byrow=T)
  EVER <- matrix(c(Fs[2], Fa[2], Sa[2], Sa[2]), Ns, Ns, byrow=T)
  KRV <-  matrix(c(Fs[3], Fa[3], Sa[3], Sa[3]), Ns, Ns, byrow=T)
  OKEE <- matrix(c(Fs[4], Fa[4], Sa[4], Sa[4]), Ns, Ns, byrow=T)
  PP <- matrix(c(Fs[5], Fa[5], Sa[5], Sa[5]), Ns, Ns, byrow=T)
  SJM <- matrix(c(Fs[6], Fa[6], Sa[6], Sa[6]), Ns, Ns, byrow=T)
  
  #combine: number of blocks = number of sites
  B <- rbind(cbind(EAST,B0,B0,B0,B0,B0),
             cbind(B0,EVER,B0,B0,B0,B0),
             cbind(B0,B0,KRV,B0,B0,B0),
             cbind(B0,B0,B0,OKEE,B0,B0),
             cbind(B0,B0,B0,B0,PP,B0),
             cbind(B0,B0,B0,B0,B0,SJM))
  return(B)
}

getB_hindcast_floater <- function(x, MUspei){
  # x = a data frame with columns for MU, age_cat, spei_cat, and demographic params
  # MUspei = dataframe with SPEI category for each MU from which to draw demographic rates
  
  
  B0 <- matrix(0, Ns, Ns)  						             # empty blocks
  
  #filter by MUspei
  tmp<- merge(MUspei, x, by = c("MU", "spei_cat"))
  
  #filter by age class
  tmp.a<- tmp %>% 
    filter(age_cat == "adult")
  tmp.j<- tmp %>% 
    filter(age_cat == "juv")
  
  #main code convert 0 to 1e-5 to get reliable values from beta (won't take 0)
  
  #define params
  Sj <- getBetaDist(meanx = tmp.j$S, stdevx = tmp.j$S - tmp.j$S_LCL)
  Sa <- getBetaDist(meanx = tmp.a$S, stdevx = tmp.a$S - tmp.a$S_LCL)
  NS <- getBetaDist(meanx = tmp.a$NS, stdevx = tmp.a$NS_UCL - tmp.a$NS)
  Jbreed <- getBetaDist(meanx = tmp.j$Breed, stdevx = tmp.j$Breed.SE)
  Abreed <- getBetaDist(meanx = tmp.a$Breed, stdevx = tmp.a$Breed.SE)
  Jattempt <- getlogNormDist(meanx = tmp.j$Attempt, stdevx = tmp.j$Attempt.SE)
  Aattempt <- getlogNormDist(meanx = tmp.a$Attempt, stdevx = tmp.a$Attempt.SE)
  Y <- getlogNormDist(meanx = tmp.a$Yng_cnt, stdevx = tmp.a$Yng_cnt.SE)
  
  #adjust very small demographic rates to be zero
  Sj[Sj < 1e-5] <- 0
  Sa[Sa < 1e-5] <- 0
  NS[NS < 1e-5] <- 0
  Jbreed[Jbreed < 1e-5] <- 0
  Abreed[Abreed < 1e-5] <- 0
  
  #fecundity: set to 0 breeding probability
  ms <- 0 * NS * Y/2 * Jattempt
  ma <- 0 * NS * Y/2 * Aattempt
  
  #fertility
  Fs <- Sj * ms
  Fa <- Sj * ma
  
  EAST <- matrix(c(Fs[1], Fa[1], Sa[1], Sa[1]), Ns, Ns, byrow=T)
  EVER <- matrix(c(Fs[2], Fa[2], Sa[2], Sa[2]), Ns, Ns, byrow=T)
  KRV <-  matrix(c(Fs[3], Fa[3], Sa[3], Sa[3]), Ns, Ns, byrow=T)
  OKEE <- matrix(c(Fs[4], Fa[4], Sa[4], Sa[4]), Ns, Ns, byrow=T)
  PP <- matrix(c(Fs[5], Fa[5], Sa[5], Sa[5]), Ns, Ns, byrow=T)
  SJM <- matrix(c(Fs[6], Fa[6], Sa[6], Sa[6]), Ns, Ns, byrow=T)
  
  #combine: number of blocks = number of sites
  B <- rbind(cbind(EAST,B0,B0,B0,B0,B0),
             cbind(B0,EVER,B0,B0,B0,B0),
             cbind(B0,B0,KRV,B0,B0,B0),
             cbind(B0,B0,B0,OKEE,B0,B0),
             cbind(B0,B0,B0,B0,PP,B0),
             cbind(B0,B0,B0,B0,B0,SJM))
  return(B)
}

#---------------------------


# Calculates M matrix using transition rates among MUs
getM <- function(disp, disp.SE){
  M0 <- matrix(0, Np, Np) # empty blocks
  
  dispmat <- as.matrix(disp)
  dispmat.SE <- as.matrix(disp.SE)
  
  #convert 0 to 0.0001 to get reliable values from beta (won't take 0)
  dispmat[dispmat < 1e-5] <- 1e-5  
  dispmat.SE[dispmat.SE < 1e-5] <- 1e-5  
  dispmat.SE[] <- 0.001
  #beta distribution
  dispersal <- getBetaDist(meanx = as.vector(dispmat), stdevx = as.vector(dispmat.SE))#vectors by column
  disp.mat <- matrix(dispersal, nrow = Np, ncol = Np, byrow=FALSE)
  
  #format dispersal matrix
  disp.mat[disp.mat < 1e-5] <- 0  #converts small values to 0
  
  #disp.mat.rsum <- rowSums(disp.mat)
  
  #should be colSums
  disp.mat.csum <- colSums(disp.mat)
  disp.mat <- disp.mat/disp.mat.csum #adjust so rows sum to 1
  
  #combine: number of blocks = number of stages
  M1 <- rbind(cbind(disp.mat,M0),
              cbind(M0,disp.mat))
  return(M1)
}




# Calculates M matrix using transition rates among MUs using two different transition matrices, one for adults and one for juveniles
getM2 <- function(disp_ad, disp_ad.SE, disp_juv, disp_juv.SE) {
  # Create an empty matrix for off-diagonal blocks (size Np x Np)
  M0 <- matrix(0, Np, Np)
  
  ## Process adult dispersal parameters:
  dispmat_ad <- as.matrix(disp_ad)
  dispmat_ad.SE <- as.matrix(disp_ad.SE)
  # Replace very small values with 1e-5 for reliability
  dispmat_ad[dispmat_ad < 1e-5] <- 1e-5  
  dispmat_ad.SE[dispmat_ad.SE < 1e-5] <- 1e-5  
  # Here we force the SE to a constant value if desired (as in your original)
  dispmat_ad.SE[] <- 0.001  
  # Compute the beta distribution draws for adult dispersal
  dispersal_ad <- getBetaDist(meanx = as.vector(dispmat_ad), stdevx = as.vector(dispmat_ad.SE))
  disp.mat.ad <- matrix(dispersal_ad, nrow = Np, ncol = Np, byrow = FALSE)
  disp.mat.ad[disp.mat.ad < 1e-5] <- 0
  # Normalize so that columns (or rows) sum to 1; here we use column sums:
  disp.mat.ad.csum <- colSums(disp.mat.ad)
  disp.mat.ad <- disp.mat.ad / disp.mat.ad.csum
  
  ## Process juvenile dispersal parameters:
  dispmat_juv <- as.matrix(disp_juv)
  dispmat_juv.SE <- as.matrix(disp_juv.SE)
  dispmat_juv[dispmat_juv < 1e-5] <- 1e-5
  dispmat_juv.SE[dispmat_juv.SE < 1e-5] <- 1e-5
  dispmat_juv.SE[] <- 0.001
  dispersal_juv <- getBetaDist(meanx = as.vector(dispmat_juv), stdevx = as.vector(dispmat_juv.SE))
  disp.mat.juv <- matrix(dispersal_juv, nrow = Np, ncol = Np, byrow = FALSE)
  disp.mat.juv[disp.mat.juv < 1e-5] <- 0
  disp.mat.juv.csum <- colSums(disp.mat.juv)
  disp.mat.juv <- disp.mat.juv / disp.mat.juv.csum
  
  # Construct the stage-structured dispersal matrix.
  # Here we assume the first block (top-left) applies to juveniles and the second (bottom-right) to adults.
  # Off-diagonals are 0 (i.e. no direct stage change during dispersal).
  M1 <- rbind(
    cbind(disp.mat.juv, M0),  # Juvenile block for juveniles, zeros for adult transition
    cbind(M0,       disp.mat.ad)  # Adult block for adults, zeros for juvenile transition
  )
  
  return(M1)
}


#---------------------------

# Wrapper for rbeta(), but takes mean and SD as input
getBetaDist <- function(meanx, stdevx){
  #adjust var if needed
  mean.check <- meanx*(1-meanx)
  for(i in 1:length(stdevx)){
    if(stdevx[i]^2 > mean.check[i]) stdevx[i] <- mean.check[i]
  }
  
  alpha <- ((1 - meanx) / stdevx^2 - 1 /meanx) * meanx ^ 2
  beta <- alpha * (1 / meanx - 1)
  alpha[is.nan(alpha)] <- 0
  beta[is.nan(beta)] <- 0
  betar <- unlist(Map(rbeta, n = 1, alpha, beta))
  return(betar)
}

#---------------------------

# Wrapper for rlnorm(), but takes meand and SD as input
getlogNormDist <- function(meanx, stdevx){
  logmean <- log(meanx) - 1/2*log((stdevx/meanx)^2 + 1)
  logsd <- (log((stdevx/meanx)^2 + 1))^0.5
  logmean[is.nan(logmean)] <- 0
  logsd[is.nan(logsd)] <- 0
  logNormr <- unlist(Map(rlnorm, n = 1, meanlog = logmean, sdlog = logsd))
  return(logNormr)
}

#-----------------------------

# Function to create flat violin plot in 'ggplot2'

# copied from
# https://gist.github.com/dgrtwo/eb7750e74997891d7c20

# somewhat hackish solution to:
# https://twitter.com/EamonCaddigan/status/646759751242620928
# based mostly on copy/pasting from ggplot2 geom_violin source:
# https://github.com/hadley/ggplot2/blob/master/R/geom-violin.r



"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(ymin = min(y),
                     ymax = max(y),
                     xmin = x,
                     xmax = x + width / 2)
          },
          
          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data, xminv = x,
                              xmaxv = x + violinwidth * (xmax - x))
            
            # Make sure it's sorted properly to draw the outline
            newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
                             plyr::arrange(transform(data, x = xmaxv), -y))
            
            # Close the polygon: set first and last point the same
            # Needed for coord_polar and such
            newdata <- rbind(newdata, newdata[1,])
            
            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },
          
          draw_key = draw_key_polygon,
          
          default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                            alpha = NA, linetype = "solid"),
          
          required_aes = c("x", "y")
  )