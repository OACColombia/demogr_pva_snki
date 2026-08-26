
# Internal functions to change the time scale of the environ covars 
#(Maru´s comment: since every variable has it unique name, I had to create one function per variable. There must  be a way to create a single function for all variables for sure)
daily2monthly_pr = function(stack, fun) {
  # stack = a RasterStack on which to convert var from daily to monthly scale (per year)
  # fun = a function specifying how to aggregate the daily data to monthly scale
  
  #get the date from the names of the layers and extract the year-month
  ind  <- sub("pr_(\\d{4}-\\d{2}).*", "\\1", names(stack))
  # ind <- as.numeric(ind)
  
  
  # Apply the function using `tapp` with the grouping vector
  tmp <- tapp(stack, ind, fun = fun)
  
  names(tmp)<- gsub("X", "pr_", names(tmp))
  
  return(tmp)
}

# Internal function to change the time scale of the environ covars
daily2monthly_tmmn = function(stack, fun) {
  # stack = a RasterStack on which to convert var from daily to monthly scale (per year)
  # fun = a function specifying how to aggregate the daily data to monthly scale
  
  #get the date from the names of the layers and extract the year-month
  ind  <- sub("tmmn_(\\d{4}-\\d{2}).*", "\\1", names(stack))
  # ind <- as.numeric(ind)
  
  
  # Apply the function using `tapp` with the grouping vector
  tmp <- tapp(stack, ind, fun = fun)
  
  names(tmp)<- gsub("X", "tmmn_", names(tmp))
  
  return(tmp)
}

# Internal function to change the time scale of the environ covars
daily2monthly_tmmx = function(stack, fun) {
  # stack = a RasterStack on which to convert var from daily to monthly scale (per year)
  # fun = a function specifying how to aggregate the daily data to monthly scale
  
  #get the date from the names of the layers and extract the year-month
  ind  <- sub("tmmx_(\\d{4}-\\d{2}).*", "\\1", names(stack))
  # ind <- as.numeric(ind)
  
  
  # Apply the function using `tapp` with the grouping vector
  tmp <- tapp(stack, ind, fun = fun)
  
  names(tmp)<- gsub("X", "tmmx_", names(tmp))
  
  return(tmp)
}







#Download gridMET data using the 'cimateR' package


params = c("precipitation_amount","daily_minimum_temperature", "daily_maximum_temperature")

gridMET_summary = function(AOI, varnames, startDate, endDate) {
  # AOI = an sf or sp object on which to extract climate data
  # varname = a vector of the model params to download from MACA
  # startDate = a character string of the start date in ISO format
  # endDate = a character string of the end date in ISO format
  
  # tictoc::tic()
  gridmet = getGridMET(AOI = AOI,
                       varname = varnames, 
                       startDate = startDate, 
                       endDate = endDate)
  # tictoc::toc()
  print(paste("Finished downloading data for", AOI$MU_Name))
  
  
  # Calculate sum total rainfall (mm) for precip and avg min/max temps
  clim<- vector('list', length = length(params))
  for (i in 1:length(params)) { #print(paste("param",i))
    ind<- grep(params[i], names(gridmet))
    
    # perform function (sum or mean) on each set of predictions over time period
    if (params[i] == 'precipitation_amount') {
      clim[[i]]<- purrr::map(gridmet[ind], ~daily2monthly_pr(., fun = sum)) 
      names(clim)[i]<- params[i]
      
    } else if (params[i] == 'daily_minimum_temperature') {
      clim[[i]]<- purrr::map(gridmet[ind], ~daily2monthly_tmmn(., fun = mean)) 
      names(clim)[i]<- params[i]
      
      
    } else {
      clim[[i]]<- purrr::map(gridmet[ind], ~daily2monthly_tmmx(., fun = mean)) 
      names(clim)[i]<- params[i]
      
    }  #close if-else conditional
  }  #close i
  
  # create vector of model names
  n.yrmonths<- nlyr(clim$precipitation_amount[[1]])
  
  date.seq<- seq.Date(as.Date(startDate, format = "%Y-%m-%d"),
                      as.Date(endDate, format = "%Y-%m-%d"),
                      by = "month")
  months1<- lubridate::month(date.seq)
  year1<- lubridate::year(date.seq)
  
  
  # summarize params for gridMET

  tmp.prcp<- t(exact_extract(rast(clim$precipitation_amount), AOI, "mean"))
  tmp.tmin<- t(exact_extract(rast(clim$daily_minimum_temperature), AOI, "mean"))
  tmp.tmax<- t(exact_extract(rast(clim$daily_maximum_temperature), AOI, "mean"))
  clim.res<- data.frame(month = months1, year = year1, prcp = tmp.prcp,
                          tmin = tmp.tmin, tmax = tmp.tmax, row.names = NULL)
    
  
  
  return(clim.res)
}


