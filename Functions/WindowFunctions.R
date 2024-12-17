
# this script contains functions for running sliding windows to identify the seasonal period over which climate variables are improtant for abundance 

library(ncdf4)
library(lubridate)
library(lme4)
library(partR2)
library(DHARMa)

options(warn=1) # print warnings as they occur

#########################################
#########################################

# Arguments for the below functions for sliding windows are as follows:
#   SpeciesName - scientific name in the form "Genus_species"
#   AbundanceData - data frame of abundance data (UKBMS site indices) 
#   FocWindow - the window you're looking for "MaxTemp1" or "Precip1"
#   LastWindow - if a previous window has been run, specify it here and it will also be included in the models
#   Scaled - TRUE/FALSE, should climate variables be scaled and centred?
#   Rows - which models to run. allows running in parallel, c(start,end) can specify "end"
#   CountPercentCurr - the day when window search will end (in the current year), given as a precentage of counts observed
#   CountPercentPrev - the day the window search will start (in the previous year), given as a precentage of counts observed (-30 means 30 days before first count)
#   WindowDurations - the durations of windows to explore
#   WindowStartInterval- the interval between window starts
#   SpAvSubset - the years spatial average should be calculated over (if not specified, all years are used)
#   AR1 - TRUE/FALSE to include previous year's abundance as a predictor
#   AR1YearDev - TRUE/FALSE - is previous year's abundance predictor using the yearly deviation from the site average abundance?
#   SpaceOnly/TimeOnly - conduct a window search for only spatial/temporal climate variation
#   StartVals - TRUE/FALSE. Do we want to use starting from the best model identified previously? THis is used when rerunning window searches where some models didn't fit well. Will write over previous run.

#########################################
#########################################

### function to get the climate variables for the windows we're exploring 

GetWindowsClimates<-function(SpeciesName, AbundanceData, CountPercentCurr, CountPercentPrev,
                             WindowDurations, WindowStartInterval, SpAvSubset=NA){

  ### Get climate values for each site in each year
 
  # get the day we're working from, when a given % of counts have been made
  SppDays<-read.csv(paste0("Data/SortedSpeciesData/",SpeciesName,"_PercentDays.csv"))
  FocDayCurr<-SppDays$Day[which(SppDays$PercentObserved==CountPercentCurr)]
  FocDayPrev<-SppDays$Day[which(SppDays$PercentObserved==CountPercentPrev)]
  
  # find the years, sites and site-years 
  Years<-unique(AbundanceData$Year)
  Sites<-unique(AbundanceData$Site)
  
  # want climate for all siteyears, including those where abundance wasn't counted (for within subject centring)
  SiteYears<-expand.grid(Sites,Years)
  colnames(SiteYears)<-c("Site","Year")
  SiteYears$siteyear<-paste(SiteYears$Site,SiteYears$Year,sep="_")
  nsiteyears<-nrow(SiteYears)
  
  # Get climates for each siteyear

  MaxTemp<-read.csv("Data/Climate/Temperature/SiteYearMaxTemp.csv")
  Precip<-read.csv("Data/Climate/Precipitation/SiteYearPrecip.csv") 
  
  # find days of climate data to consider in the current and previous years
  # this is working backwards, starting with FocDayCurr as day 1
  NCurrYearDays<-length(1:FocDayCurr) 
  PrevYearStart<-FocDayPrev+1 # day of previous year to start from - this is the day after CountPercentPrev % of individuals have been observed
  NPrevYearDays<-length(PrevYearStart:365) 
  # because we're looking at the days prior to FocDayCurr, prev year ordinal day numbers will be shifted if prev year is a leap year
  
  ndays<-NPrevYearDays+NCurrYearDays 
  
  DayNumbers<-data.frame(DaysPrior=ndays:1,OrdinalDay=c(PrevYearStart:365,1:FocDayCurr),
                         OrdinalDayLeapPrevYear=c((PrevYearStart+1):366,1:FocDayCurr),
                         Year=c(rep("Previous",NPrevYearDays),rep("Current",NCurrYearDays)),
                         Date=rep(NA,ndays))
  
  DayNumbers$Date<-paste0(as.POSIXlt(as.character(DayNumbers$OrdinalDay),format="%j")$mday,
                          month.name[as.POSIXlt(as.character(DayNumbers$OrdinalDay),format="%j")$mon+1]) # +1 as this is months after the first month of the year 
  # dates will be the same in leap and non-leap years, but ordinal day numbers are +1 in leap year
  
  
  # Generate a covariate matrix for each climate variable to store daily climate values for each siteyear
  # This has a row for each siteyear and a column for each day of climate we want to consider
  # day numbers are days prior to FocDay + 1, i.e. FocDay is day 1
  CovariateMatrixMaxTemp<-CovariateMatrixPrecip<-matrix(nrow=nsiteyears,ncol=ndays)
  rownames(CovariateMatrixMaxTemp)<-rownames(CovariateMatrixPrecip)<-SiteYears$siteyear
  colnames(CovariateMatrixMaxTemp)<-colnames(CovariateMatrixPrecip)<-ndays:1
  
  # fill out covariate matrices with climate data
  
  for(focrow in 1:nrow(CovariateMatrixMaxTemp)){
    
    # find the siteyear, site and year for this row
    focsiteyear<-rownames(CovariateMatrixMaxTemp)[focrow] 
    focsite<-as.numeric(strsplit(focsiteyear,split="_")[[1]][1])
    focyear<-as.numeric(strsplit(focsiteyear,split="_")[[1]][2])
    
    print(focsiteyear)
    
    # find previous year and siteyear
    prevyear<-focyear-1
    prevsiteyear<-paste0(focsite,"_",prevyear)
    
    # get current and previous year's climate  
    
    # the days used depends on leap years
    
    # day numbers of interest are based on when a proportion of counts have been observed across all sites and years
    # this gives the ordinal day number rather than a specific date, so it doesn't matter if the current year is a leap year
    # and we don't need to +1 to the FocDayCurr (which we would if we were interested in a specific date)
    # we only have to do something if the previous year is a leap year because we are considering days prior to FocDayCurr 
    # so this will change which days' climate we take for the previous year
    
    if(leap_year(prevyear)==FALSE){ # if previous year is not a leap year
      
      # first NPrevYearDays columns of covariate matrices are filled by previous year's climate
      # climate values are for days PrevYearStart:365 of previous year
      CovariateMatrixMaxTemp[focrow,c(1:NPrevYearDays)]<-as.numeric(MaxTemp[which(MaxTemp$SiteYear==prevsiteyear),
                                                                            match(paste0("Day",PrevYearStart:365),colnames(MaxTemp))])
      
      CovariateMatrixPrecip[focrow,c(1:NPrevYearDays)]<-as.numeric(Precip[which(Precip$SiteYear==prevsiteyear),
                                                                          match(paste0("Day",PrevYearStart:365),colnames(Precip))])
      
      # the rest of the columns, from NPrevYearDays+1 to ndays, are the current year
      # climate values are for days 1:FocDayCurr  
      CovariateMatrixMaxTemp[focrow,c((NPrevYearDays+1):ndays)]<-as.numeric(MaxTemp[which(MaxTemp$SiteYear==focsiteyear),
                                                                                    match(paste0("Day",1:FocDayCurr),colnames(MaxTemp))])
      
     CovariateMatrixPrecip[focrow,c((NPrevYearDays+1):ndays)]<-as.numeric(Precip[which(Precip$SiteYear==focsiteyear),
                                                                                  match(paste0("Day",1:FocDayCurr),colnames(Precip))])
      
    }else{ # if previous year is a leap year
      
      # this changes which values we take for the previous year
      
      # leap year in previous year shifts ordinal dates by a day
      # so last year's climate date starts one day later at PrevYearStart+1 and ends on day 366
      # this is because there's an extra day in the previous year
      # and we're only looking at days prior to FocDayCurr
      CovariateMatrixMaxTemp[focrow,c(1:NPrevYearDays)]<-as.numeric(MaxTemp[which(MaxTemp$SiteYear==prevsiteyear),
                                                                            match(paste0("Day",(PrevYearStart+1):366),colnames(MaxTemp))])

      CovariateMatrixPrecip[focrow,c(1:NPrevYearDays)]<-as.numeric(Precip[which(Precip$SiteYear==prevsiteyear),
                                                                          match(paste0("Day",(PrevYearStart+1):366),colnames(Precip))])
      
      # the rest of the columns, from NPrevYearDays+1 to ndays, are the current year
      # climate values are for days 1:FocDayCurr of the current year 
      CovariateMatrixMaxTemp[focrow,c((NPrevYearDays+1):ndays)]<-as.numeric(MaxTemp[which(MaxTemp$SiteYear==focsiteyear),
                                                                                    match(paste0("Day",1:FocDayCurr),colnames(MaxTemp))])
      
      CovariateMatrixPrecip[focrow,c((NPrevYearDays+1):ndays)]<-as.numeric(Precip[which(Precip$SiteYear==focsiteyear),
                                                                                  match(paste0("Day",1:FocDayCurr),colnames(Precip))])
      
    }
  }

  ### Find the windows we want to explore
  
  # Remove siteyears with missing climate data
  siteyearsdrop<-unique(c(rownames(which(is.na(CovariateMatrixMaxTemp),arr.ind=TRUE)),
                          rownames(which(is.na(CovariateMatrixPrecip),arr.ind=TRUE))))

  if(length(siteyearsdrop)>0){
  AbundanceData<-AbundanceData[-match(intersect(siteyearsdrop,AbundanceData$SiteYear),AbundanceData$SiteYear),]
  }
  
  # find remaining site years combinations (including those without counts for WSC)
  FilteredSiteYears<-expand.grid(unique(AbundanceData$Site),unique(AbundanceData$Year))
  colnames(FilteredSiteYears)<-c("Site","Year")
  FilteredSiteYears$SiteYear<-paste(FilteredSiteYears$Site,FilteredSiteYears$Year, sep="_")
  nFilteredSiteYears<-nrow(FilteredSiteYears)
  
  # windows are going backwards from FocDay
  # numbers are days prior to FocDay+1 (i.e. FocDay=1), mas in DayNumbers$DaysPrior
  WindowStartDays<-c(seq(1,nrow(DayNumbers),WindowStartInterval)) # intervals for starting windows

  ClimWindows<-expand.grid(WindowDurations,WindowStartDays) # data frame of all combinations of start and duration
  colnames(ClimWindows)<-c("WindowDuration","WindowStartDay")
  ClimWindows$WindowEndDay<-ClimWindows$WindowStartDay+(ClimWindows$WindowDuration-1) # window end days
  ClimWindows<-ClimWindows[-which(ClimWindows$WindowEndDay>nrow(DayNumbers)),] # only looking back to FocDay+1 in prev year 
  ClimWindows$window<-paste0(ClimWindows$WindowStartDay,"-",ClimWindows$WindowEndDay)
  
  # add dates and ordinal days corresponding to these windows
  ClimWindows$StartDate<-DayNumbers$Date[match(ClimWindows$WindowStartDay,DayNumbers$DaysPrior)]
  ClimWindows$EndDate<-DayNumbers$Date[match(ClimWindows$WindowEndDay,DayNumbers$DaysPrior)]
  
  ClimWindows$StartOrdinalDay<-DayNumbers$OrdinalDay[match(ClimWindows$WindowStartDay,DayNumbers$DaysPrior)]
  ClimWindows$EndOrdinalDay<-DayNumbers$OrdinalDay[match(ClimWindows$WindowEndDay,DayNumbers$DaysPrior)]
  
  ClimWindows$StartOrdinalDayLeapPrevYear<-DayNumbers$OrdinalDayLeapPrevYear[match(ClimWindows$WindowStartDay,DayNumbers$DaysPrior)]
  ClimWindows$EndOrdinalDayLeapPrevYear<-DayNumbers$OrdinalDayLeapPrevYear[match(ClimWindows$WindowEndDay,DayNumbers$DaysPrior)]
  
  ClimWindows$Year<-DayNumbers$Year[match(ClimWindows$WindowEndDay,DayNumbers$DaysPrior)]
  
  # Generate data frame to store windows' mean climates for each siteyear and the scaled data
  # col for each window, row for each siteyear (this is with all site years, for WSC)
  WindowsMaxTemp<-WindowsPrecip<-as.data.frame(matrix(nrow=nrow(CovariateMatrixMaxTemp),ncol=(nrow(ClimWindows)*3)))
  WindowsMaxTempScaled<-WindowsPrecipScaled<-as.data.frame(matrix(nrow=nrow(CovariateMatrixMaxTemp),ncol=(nrow(ClimWindows)*3)))
  
  
  # name windows
  colnames(WindowsMaxTemp)<-c(paste0("MaxTemp_",ClimWindows$window),paste0("MaxTemp_",ClimWindows$window,"SpatialAv"),paste0("MaxTemp_",ClimWindows$window,"YearDev"))
  colnames(WindowsPrecip)<-c(paste0("Precip_",ClimWindows$window),paste0("Precip_",ClimWindows$window,"SpatialAv"),paste0("Precip_",ClimWindows$window,"YearDev"))
  
  colnames(WindowsMaxTempScaled)<-c(paste0("MaxTempScaled_",ClimWindows$window),paste0("MaxTempScaled_",ClimWindows$window,"SpatialAv"),paste0("MaxTempScaled_",ClimWindows$window,"YearDev"))
  colnames(WindowsPrecipScaled)<-c(paste0("PrecipScaled_",ClimWindows$window),paste0("PrecipScaled_",ClimWindows$window,"SpatialAv"),paste0("PrecipScaled_",ClimWindows$window,"YearDev"))
  
  # add site/year columns
  WindowsMaxTemp$SiteYear<-WindowsPrecip$SiteYear<-rownames(CovariateMatrixMaxTemp) 
  WindowsMaxTempScaled$SiteYear<-WindowsPrecipScaled$SiteYear<-rownames(CovariateMatrixMaxTemp) 
  
  WindowsMaxTemp$Site<-WindowsPrecip$Site<-FilteredSiteYears$Site[match(WindowsMaxTemp$SiteYear,FilteredSiteYears$SiteYear)]
  WindowsMaxTempScaled$Site<-WindowsPrecipScaled$Site<-FilteredSiteYears$Site[match(WindowsMaxTempScaled$SiteYear,FilteredSiteYears$SiteYear)]
  
  WindowsMaxTemp$Year<-WindowsPrecip$Year<-FilteredSiteYears$Year[match(WindowsMaxTemp$SiteYear,FilteredSiteYears$SiteYear)]
  WindowsMaxTempScaled$Year<-WindowsPrecipScaled$Year<-FilteredSiteYears$Year[match(WindowsMaxTempScaled$SiteYear,FilteredSiteYears$SiteYear)]
  
  # covar matrices contain daily climate values for days prior to FocDayCurr+1 (i.e. FocDayCurr is day 1), as in DayNumbers$DaysPrior
  # so just need to average these for each window. 
  
  for(WindowRow in 1:nrow(ClimWindows)){ # for each window
    
    # get window name, start and end 
    focwindow<-ClimWindows$window[WindowRow]
    focstart<-ClimWindows$WindowStartDay[WindowRow]
    focend<-ClimWindows$WindowEndDay[WindowRow]
    
    # find the mean climate values for this window and store 
    WindowsMaxTemp[,which(colnames(WindowsMaxTemp)==paste0("MaxTemp_",focwindow))]<-rowMeans(CovariateMatrixMaxTemp[,match(c(focstart:focend),colnames(CovariateMatrixMaxTemp))])
    
    WindowsPrecip[,which(colnames(WindowsPrecip)==paste0("Precip_",focwindow))]<-rowMeans(CovariateMatrixPrecip[,match(c(focstart:focend),colnames(CovariateMatrixPrecip))])
    
  }
  
  # Scale data by sd and mean centre
  
  # store for centre and scale used
  ClimWindows$MaxTempScale<-NA
  ClimWindows$MaxTempCentre<-NA
  
  ClimWindows$PrecipScale<-NA
  ClimWindows$PrecipCentre<-NA
  
  for(WindowRow in 1:nrow(ClimWindows)){
    
    focwindow<-ClimWindows$window[WindowRow]
    focstart<-ClimWindows$WindowStartDay[WindowRow]
    focend<-ClimWindows$WindowEndDay[WindowRow]
    
    #scale
    MaxTempScale<-scale(WindowsMaxTemp[,which(colnames(WindowsMaxTemp)==paste0("MaxTemp_",focwindow))])
    PrecipScale<-scale(WindowsPrecip[,which(colnames(WindowsPrecip)==paste0("Precip_",focwindow))])
    
    WindowsMaxTempScaled[,which(colnames(WindowsMaxTempScaled)==paste0("MaxTempScaled_",focwindow))]<-as.numeric(MaxTempScale)
    WindowsPrecipScaled[,which(colnames(WindowsPrecipScaled)==paste0("PrecipScaled_",focwindow))]<-as.numeric(PrecipScale)
    
    
    # store centre and scale used for each variable
    ClimWindows$MaxTempScale[WindowRow]<-attributes(MaxTempScale)$`scaled:scale`
    ClimWindows$MaxTempCentre[WindowRow]<-attributes(MaxTempScale)$`scaled:center`
    ClimWindows$PrecipScale[WindowRow]<-attributes(PrecipScale)$`scaled:scale`
    ClimWindows$PrecipCentre[WindowRow]<-attributes(PrecipScale)$`scaled:center`
  }
  
  # get vars for within subject centring
  for(i in 1:nrow(ClimWindows)){
    focwindow<-ClimWindows$window[i]
    
    # find SpatialClim, the mean for each site across years 

    if(is.na(SpAvSubset[1])){
    WindowsMaxTemp[,which(colnames(WindowsMaxTemp)==paste0("MaxTemp_",focwindow,"SpatialAv"))]<-tapply(WindowsMaxTemp[,which(colnames(WindowsMaxTemp)==paste0("MaxTemp_",focwindow))],
                      WindowsMaxTemp$Site,mean)[as.character(WindowsMaxTemp$Site)]
    WindowsPrecip[,which(colnames(WindowsPrecip)==paste0("Precip_",focwindow,"SpatialAv"))]<-tapply(WindowsPrecip[,which(colnames(WindowsPrecip)==paste0("Precip_",focwindow))],
                      WindowsPrecip$Site,mean)[as.character(WindowsPrecip$Site)]

    WindowsMaxTempScaled[,which(colnames(WindowsMaxTempScaled)==paste0("MaxTempScaled_",focwindow,"SpatialAv"))]<-tapply(WindowsMaxTempScaled[,which(colnames(WindowsMaxTempScaled)==paste0("MaxTempScaled_",focwindow))],
                      WindowsMaxTempScaled$Site,mean)[as.character(WindowsMaxTempScaled$Site)]
    WindowsPrecipScaled[,which(colnames(WindowsPrecipScaled)==paste0("PrecipScaled_",focwindow,"SpatialAv"))]<-tapply(WindowsPrecipScaled[,which(colnames(WindowsPrecipScaled)==paste0("PrecipScaled_",focwindow))],
                      WindowsPrecipScaled$Site,mean)[as.character(WindowsPrecipScaled$Site)]
    
    }else{
      WindowsMaxTemp[,which(colnames(WindowsMaxTemp)==paste0("MaxTemp_",focwindow,"SpatialAv"))]<-tapply(WindowsMaxTemp[which(is.na(match(WindowsMaxTemp$Year,SpAvSubset))==FALSE),which(colnames(WindowsMaxTemp)==paste0("MaxTemp_",focwindow))],
                                                                                                         WindowsMaxTemp$Site[which(is.na(match(WindowsMaxTemp$Year,SpAvSubset))==FALSE)],mean)[as.character(WindowsMaxTemp$Site)]
      WindowsPrecip[,which(colnames(WindowsPrecip)==paste0("Precip_",focwindow,"SpatialAv"))]<-tapply(WindowsPrecip[which(is.na(match(WindowsPrecip$Year,SpAvSubset))==FALSE),which(colnames(WindowsPrecip)==paste0("Precip_",focwindow))],
                                                                                                      WindowsPrecip$Site[which(is.na(match(WindowsPrecip$Year,SpAvSubset))==FALSE)],mean)[as.character(WindowsPrecip$Site)]
      
      WindowsMaxTempScaled[,which(colnames(WindowsMaxTempScaled)==paste0("MaxTempScaled_",focwindow,"SpatialAv"))]<-tapply(WindowsMaxTempScaled[which(is.na(match(WindowsMaxTempScaled$Year,SpAvSubset))==FALSE),which(colnames(WindowsMaxTempScaled)==paste0("MaxTempScaled_",focwindow))],
                                                                                                     WindowsMaxTempScaled$Site[which(is.na(match(WindowsMaxTempScaled$Year,SpAvSubset))==FALSE)],mean)[as.character(WindowsMaxTempScaled$Site)]
      WindowsPrecipScaled[,which(colnames(WindowsPrecipScaled)==paste0("PrecipScaled_",focwindow,"SpatialAv"))]<-tapply(WindowsPrecipScaled[which(is.na(match(WindowsPrecipScaled$Year,SpAvSubset))==FALSE),which(colnames(WindowsPrecipScaled)==paste0("PrecipScaled_",focwindow))],
                                                                                                      WindowsPrecipScaled$Site[which(is.na(match(WindowsPrecipScaled$Year,SpAvSubset))==FALSE)],mean)[as.character(WindowsPrecipScaled$Site)]
    }
    
    # find YearDev, each year's deviation from this spatial mean 
    WindowsMaxTemp[,which(colnames(WindowsMaxTemp)==paste0("MaxTemp_",focwindow,"YearDev"))]<-WindowsMaxTemp[,which(colnames(WindowsMaxTemp)==paste0("MaxTemp_",focwindow))]-
                     WindowsMaxTemp[,which(colnames(WindowsMaxTemp)==paste0("MaxTemp_",focwindow,"SpatialAv"))]
    WindowsPrecip[,which(colnames(WindowsPrecip)==paste0("Precip_",focwindow,"YearDev"))]<-WindowsPrecip[,which(colnames(WindowsPrecip)==paste0("Precip_",focwindow))]-
                      WindowsPrecip[,which(colnames(WindowsPrecip)==paste0("Precip_",focwindow,"SpatialAv"))]
    
    WindowsMaxTempScaled[,which(colnames(WindowsMaxTempScaled)==paste0("MaxTempScaled_",focwindow,"YearDev"))]<-WindowsMaxTempScaled[,which(colnames(WindowsMaxTempScaled)==paste0("MaxTempScaled_",focwindow))]-
                     WindowsMaxTempScaled[,which(colnames(WindowsMaxTempScaled)==paste0("MaxTempScaled_",focwindow,"SpatialAv"))]
    WindowsPrecipScaled[,which(colnames(WindowsPrecipScaled)==paste0("PrecipScaled_",focwindow,"YearDev"))]<-WindowsPrecipScaled[,which(colnames(WindowsPrecipScaled)==paste0("PrecipScaled_",focwindow))]-
                      WindowsPrecipScaled[,which(colnames(WindowsPrecipScaled)==paste0("PrecipScaled_",focwindow,"SpatialAv"))]
  
    }
  
  # only keep those siteyears with counts, not all those used for WSC
  WindowsMaxTemp<-WindowsMaxTemp[match(FilteredSiteYears$SiteYear,WindowsMaxTemp$SiteYear),]
  WindowsPrecip<-WindowsPrecip[match(FilteredSiteYears$SiteYear,WindowsPrecip$SiteYear),]
  WindowsMaxTempScaled<-WindowsMaxTempScaled[match(FilteredSiteYears$SiteYear,WindowsMaxTempScaled$SiteYear),]
  WindowsPrecipScaled<-WindowsPrecipScaled[match(FilteredSiteYears$SiteYear,WindowsPrecipScaled$SiteYear),]
  
  # Merge window climate values into site indices data frame
  AbundanceData<-merge(AbundanceData,WindowsMaxTemp,by=c("SiteYear","Site","Year"))
  AbundanceData<-merge(AbundanceData,WindowsPrecip,by=c("SiteYear","Site","Year"))
  AbundanceData<-merge(AbundanceData,WindowsMaxTempScaled,by=c("SiteYear","Site","Year"))
  AbundanceData<-merge(AbundanceData,WindowsPrecipScaled,by=c("SiteYear","Site","Year"))
  

    save(AbundanceData, ClimWindows, DayNumbers, CovariateMatrixMaxTemp, CovariateMatrixPrecip, SpAvSubset,
         file=paste0("SlidingWindows/",SpeciesName,"/", CountPercentPrev,"-",CountPercentCurr,"/",SpeciesName,"_",
                     "Windows_",WindowStartInterval,"Int_",
                     WindowDurations[1],"-",WindowDurations[length(WindowDurations)],"Dur",
                     if(!is.na(SpAvSubset[1])){paste0("_",paste0(SpAvSubset[1],"-",SpAvSubset[length(SpAvSubset)]),"SpAv")},".rda"))
}



#########################################
#########################################

# function to run sliding window searches 

RunSlidingWindow<-function(SpeciesName, FocWindow, LastWindow=FALSE, 
                           Scaled=FALSE, Rows="All", CountPercentCurr, CountPercentPrev,
                           WindowDurations, WindowStartInterval, SpAvSubset=NA, AR1=FALSE, AR1YearDev=FALSE,
                           SpaceOnly=FALSE,TimeOnly=FALSE, StartVals=FALSE){
 
  # specify folder and file name for saving, dependent on function arguments
  modfolder<-paste0("SlidingWindows/",SpeciesName,"/",CountPercentPrev,"-",CountPercentCurr,"/",
                    if(AR1==TRUE){paste0("AR1/")},
                    if(SpaceOnly==TRUE){paste0("Space/")},
                    if(TimeOnly==TRUE){paste0("Time/")})
  
  modname<-paste0(if(Scaled==TRUE){paste0("_Scaled")}, # scaled?
                  "_",WindowStartInterval,"Int_",WindowDurations[1],"-",WindowDurations[length(WindowDurations)],"Dur", # start interval and duration
                  if(AR1==TRUE){paste0("_AR1")}, # AR1?
                  if(AR1YearDev==TRUE){paste0("YearDev")},
                  if(!is.na(SpAvSubset[1])){paste0("_",paste0(SpAvSubset[1],"-",SpAvSubset[length(SpAvSubset)]),"SpAv")}, # SpatialAv subset?
                  if(SpaceOnly==TRUE){paste0("_Space")}, if(TimeOnly==TRUE){paste0("_Time")})# space or time alone?
                  
  
  # find the day we're working from, when  of counts have been made
  SppDays<-read.csv(paste0("Data/SortedSpeciesData/",SpeciesName,"_PercentDays.csv"))

  # load windows to search
  load(paste0("SlidingWindows/",SpeciesName,"/", CountPercentPrev,"-",CountPercentCurr,"/",SpeciesName,"_",
              "Windows_",WindowStartInterval,"Int_",
              WindowDurations[1],"-",WindowDurations[length(WindowDurations)],"Dur",
              if(!is.na(SpAvSubset[1])){paste0("_",paste0(SpAvSubset[1],"-",SpAvSubset[length(SpAvSubset)]),"SpAv")},".rda"))

  
  # specify if we are inputting scaled variables or not
  if(Scaled==TRUE){
    MaxTempVar<-"MaxTempScaled"
    PrecipVar<-"PrecipScaled"
  }else{
    MaxTempVar<-"MaxTemp"
    PrecipVar<-"Precip"
  }
  
  # is focal window temperature or precipitation?
  # needed for extracting model estimates and setting start values
  if(length(grep("MaxTemp",FocWindow))==1){focVar<-MaxTempVar}
  if(length(grep("Precip",FocWindow))==1){focVar<-PrecipVar}
  
  # if AR1 get previous year's abundances
  if(AR1==TRUE){
  AbundanceData$PrevYrSITE.INDEX<-NA
  for(focrow in 1:nrow(AbundanceData)){
    focsite<-AbundanceData$Site[focrow]
    prevyear<-as.numeric(paste0(AbundanceData$Year[focrow]))-1
    prevrow<-intersect(which(AbundanceData$Site==focsite),which(AbundanceData$Year==prevyear))
    if(length(prevrow)>0){
      AbundanceData$PrevYrSITE.INDEX[focrow]<-AbundanceData$SITE.INDEX[prevrow]
    }
  }
  AbundanceData<-AbundanceData[-which(is.na(AbundanceData$PrevYrSITE.INDEX)),] # remove those where previous year's count is not available
  AbundanceData$logPrevYrSITE.INDEX<-log(AbundanceData$PrevYrSITE.INDEX+1)
  
  if(AR1YearDev==TRUE){
    AbundanceData$logPrevYrSITE.INDEXSpatialAv<-tapply(AbundanceData$logPrevYrSITE.INDEX,AbundanceData$Site,mean)[as.character(AbundanceData$Site)]
    AbundanceData$logPrevYrSITE.INDEXYearDev<-AbundanceData$logPrevYrSITE.INDEX-AbundanceData$logPrevYrSITE.INDEXSpatialAv
   }
  }
  
  print(paste0("n = ",nrow(AbundanceData)))

  if(LastWindow==FALSE){ # if this is the first window set up store for best windows
    
    BestWindows<-data.frame(matrix(nrow=3,ncol=2))
    rownames(BestWindows)<-c("Window","AIC","AICChange")
    colnames(BestWindows)<-c("MaxTemp1","Precip1") 
    PrevWindows<-c()
    
  }else{  # if this is not the first window, load the previous window(s)
    
      load(paste0(modfolder,SpeciesName,"_SlidingWindowOutputs",modname,if(StartVals==TRUE){paste0("_StartVals")},".rda"))
      
    # which windows were previously run?
    PrevWindows<-colnames(BestWindows)[which(is.na(BestWindows["Window",])==FALSE)]
  }
  
  
  WindowModels<-data.frame(windows=ClimWindows$window) # data frame of possible temp windows
  colnames(WindowModels)[1]<-paste0("Window",FocWindow)
  
  # add any previously identified windows to table
  if(length(intersect(PrevWindows,"MaxTemp1"))==1){WindowModels$WindowMaxTemp1<-BestWindows["Window","MaxTemp1"]}
  if(length(intersect(PrevWindows,"Precip1"))==1){WindowModels$WindowPrecip1<-BestWindows["Window","Precip1"]}

  # Second window cannot overlap with first window for a variable
  # so remove potential windows that overlap with window 1

  nwindows<-nrow(WindowModels)
  
  WindowModels$AIC<-NA
  WindowModels$logLik<-NA
  
  if(SpeciesName=="Anthocharis_cardamines"){
    WindowModels$R2M<-NA
    WindowModels$R2C<-NA
  }
  
  WindowModels$InterceptEst<-NA
  WindowModels$InterceptSE<-NA
  
  if(AR1==TRUE){
    if(AR1YearDev==TRUE){
      WindowModels$logPrevYrSITE.INDEXYearDevEst<-NA
      WindowModels$logPrevYrSITE.INDEXYearDevSE<-NA
    }else{
  WindowModels$logPrevYrSITE.INDEXEst<-NA
  WindowModels$logPrevYrSITE.INDEXSE<-NA
    }
  }
  
  # add columns for the focal window estimates
  ncol<-ncol(WindowModels)
  
  if(SpaceOnly==FALSE&&TimeOnly==FALSE){
  WindowModels[,(ncol+1):(ncol+8)]<-NA
  colnames(WindowModels)[(ncol+1):(ncol+8)]<-c(paste0(FocWindow, c("SpatialAvEst","SpatialAvSE",
          "SpatialAvSqEst","SpatialAvSqSE","YearDevEst","YearDevSE","SpatialAvYearDevEst",
          "SpatialAvYearDevSE")))
  
  # if previous windows have run, add columns for these and their estimates
  if(LastWindow!=FALSE){
    ncol<-ncol(WindowModels)
    WindowModels[,((ncol+1):(ncol+(length(PrevWindows)*8)))]<-NA
    colnames(WindowModels)[((ncol+1):(ncol+(length(PrevWindows)*8)))]<-as.vector(outer(PrevWindows, c("SpatialAvEst","SpatialAvSE",
        "SpatialAvSqEst","SpatialAvSqSE","YearDevEst","YearDevSE","SpatialAvYearDevEst",
        "SpatialAvYearDevSE"),paste0))
    }
  }
  
  if(SpaceOnly==TRUE){
    WindowModels[,(ncol+1):(ncol+4)]<-NA
    colnames(WindowModels)[(ncol+1):(ncol+4)]<-c(paste0(FocWindow, c("SpatialAvEst","SpatialAvSE","SpatialAvSqEst","SpatialAvSqSE"))) 
    
    # if previous windows have run, add columns for these and their estimates
    if(LastWindow!=FALSE){
      ncol<-ncol(WindowModels)
      WindowModels[,((ncol+1):(ncol+(length(PrevWindows)*4)))]<-NA
      colnames(WindowModels)[((ncol+1):(ncol+(length(PrevWindows)*4)))]<-as.vector(outer(PrevWindows, c("SpatialAvEst","SpatialAvSE","SpatialAvSqEst","SpatialAvSqSE"),paste0))
    }
  }
  
  if(TimeOnly==TRUE){
    WindowModels[,(ncol+1):(ncol+4)]<-NA
    colnames(WindowModels)[(ncol+1):(ncol+4)]<-c(paste0(FocWindow, c("YearDevEst","YearDevSE", 
                                                                     "SpatialAvYearDevEst","SpatialAvYearDevSE")))
    # if previous windows have run, add columns for these and their estimates
    if(LastWindow!=FALSE){
      ncol<-ncol(WindowModels)
      WindowModels[,((ncol+1):(ncol+(length(PrevWindows)*4)))]<-NA
      colnames(WindowModels)[((ncol+1):(ncol+(length(PrevWindows)*4)))]<-as.vector(outer(PrevWindows, c(
        "YearDevEst","YearDevSE","SpatialAvYearDevEst","SpatialAvYearDevSE"),paste0))
    }
  }

    WindowModels$ConvergeWarning<-WindowModels$SingularWarning<-WindowModels$EigenvalueWarning<-NA
  
  # find the columns for the windows to run
  tempcols<-grep("WindowMaxTemp",colnames(WindowModels))
  precipcols<-grep("WindowPrecip",colnames(WindowModels))
  
  # if we're re-running with start values, load the previous window run and get the estimates for the best model
  if(StartVals==TRUE){
    
    # load previous run
    PrevWindowModels<-read.table(paste0(modfolder,SpeciesName,"_SlidingWindowRun",FocWindow,modname,".txt"),header = TRUE)
    
    # columns for the best windows
    besttempcols<-grep("WindowMaxTemp",colnames(PrevWindowModels))
    bestprecipcols<-grep("WindowPrecip",colnames(PrevWindowModels))
    
    # row for best model
    bestrow<-which(PrevWindowModels$AIC==min(PrevWindowModels$AIC,na.rm=T))
    
      if(SpaceOnly==FALSE&&TimeOnly==FALSE){
        bestformula<-as.formula(paste0("SITE.INDEX~",if(AR1==TRUE){if(AR1YearDev==TRUE){paste0("logPrevYrSITE.INDEXYearDev+")}else{paste0("logPrevYrSITE.INDEX+")}},
                                      if(length(besttempcols)>0){paste0("`",MaxTempVar,"_",PrevWindowModels[bestrow,besttempcols],"SpatialAv`*`",MaxTempVar,"_",PrevWindowModels[bestrow,besttempcols],"YearDev`+I(`",MaxTempVar,"_",PrevWindowModels[bestrow,besttempcols],"SpatialAv`^2)",collapse="",sep="+")},                        
                                      if(length(bestprecipcols)>0){paste0("`",PrecipVar,"_",PrevWindowModels[bestrow, bestprecipcols], "SpatialAv`*`",PrecipVar,"_",PrevWindowModels[bestrow, bestprecipcols],"YearDev`+I(`",PrecipVar,"_",PrevWindowModels[bestrow, bestprecipcols],"SpatialAv`^2)",sep = "+",collapse="")},
                                      "(1|Year)+(1|Site)+(1|GridCell50km)+(1|GridCell50km:Year)+(1|resid)"))
      }
    
    if(SpaceOnly==TRUE){
      bestformula<-as.formula(paste0("SITE.INDEX~",if(AR1==TRUE){if(AR1YearDev==TRUE){paste0("logPrevYrSITE.INDEXYearDev+")}else{paste0("logPrevYrSITE.INDEX+")}},
                                    if(length(besttempcols)>0){paste0("`",MaxTempVar,"_",PrevWindowModels[bestrow,besttempcols],"SpatialAv`+I(`",MaxTempVar,"_",PrevWindowModels[bestrow,besttempcols],"SpatialAv`^2)",collapse="",sep="+")},                        
                                    if(length(bestprecipcols)>0){paste0("`",PrecipVar,"_",PrevWindowModels[bestrow, bestprecipcols], "SpatialAv`+I(`",PrecipVar,"_",PrevWindowModels[bestrow, bestprecipcols],"SpatialAv`^2)",sep = "+",collapse="")},
                                    "(1|Year)+(1|Site)+(1|GridCell50km)+(1|GridCell50km:Year)+(1|resid)"))
    }
    
    if(TimeOnly==TRUE){
      bestformula<-as.formula(paste0("SITE.INDEX~",if(AR1==TRUE){if(AR1YearDev==TRUE){paste0("logPrevYrSITE.INDEXYearDev+")}else{paste0("logPrevYrSITE.INDEX+")}},
                                    if(length(besttempcols)>0){paste0("`",MaxTempVar,"_",PrevWindowModels[bestrow,besttempcols],"YearDev`+`",MaxTempVar,"_",PrevWindowModels[bestrow,besttempcols],"YearDev`:`",MaxTempVar,"_",PrevWindowModels[bestrow,besttempcols],"SpatialAv`",collapse="",sep="+")},                        
                                    if(length(bestprecipcols)>0){paste0("`",PrecipVar,"_",PrevWindowModels[bestrow, bestprecipcols],"YearDev`+`",PrecipVar,"_",PrevWindowModels[bestrow, bestprecipcols], "YearDev`:`",PrecipVar,"_",PrevWindowModels[bestrow, bestprecipcols],"SpatialAv`",sep = "+",collapse="")},
                                    "(1|Year)+(1|Site)+(1|GridCell50km)+(1|GridCell50km:Year)+(1|resid)"))
    }
     
       # if using Site Indices, we want a Poisson model
      bestmod<-try(glmer(bestformula,family=poisson,data=AbundanceData,na.action=na.omit))
      
      if(class(bestmod)=="try-error"){
                bestmod<-glmer(bestformula,family=poisson,data=AbundanceData,na.action=na.omit,glmerControl(optimizer="bobyqa"))
      }
    
    # get start values - these are the estimates from the best window of the last model run
    StartValues<- fixef(bestmod)
  }
  
  # find rows to run
  # this allows different rows to be run in parallel
  if(any(Rows=="All")){
    RowsToRun<-1:nrow(WindowModels)
  }else{
    if((any(Rows=="end"))){
      RowsToRun<-Rows[1]:nrow(WindowModels)
    }else{
      if(length(Rows)>1){
      RowsToRun<-Rows[1]:Rows[2]
      }else{
        RowsToRun<-Rows
      }
    }
  }
  
  # load any models that have already been run
  if(file.exists(paste0(modfolder,SpeciesName,"_SlidingWindowRun",FocWindow,modname,if(StartVals==TRUE){paste0("_StartVals")},".txt"))){
  WindowModelsRUN<-read.table(paste0(modfolder,SpeciesName,"_SlidingWindowRun",FocWindow,modname,if(StartVals==TRUE){paste0("_StartVals")},".txt"),header = TRUE)
  }else{
    WindowModelsRUN<-NA
  }
  
  # fit models
  for(focrow in RowsToRun){
    
    print(paste0(focrow,"/",max(RowsToRun)))
    
    # if this window hasn't already run, run it 
    if(!file.exists(paste0(modfolder,SpeciesName,"_SlidingWindowRun",FocWindow,modname,if(StartVals==TRUE){paste0("_StartVals")},".txt"))|| # if there is no file yet, or...
      if(any(!is.na(WindowModelsRUN))){length(which(WindowModelsRUN[,paste0("Window",FocWindow)]==WindowModels[focrow,paste0("Window",FocWindow)]))==0}else{TRUE}){
    # ...if this window is not in the file (else clause is for if all the entries are NA)
      
    # specify model formula from table and run model
    
      if(SpaceOnly==FALSE&&TimeOnly==FALSE){
      focformula<-as.formula(paste0("SITE.INDEX~",if(AR1==TRUE){if(AR1YearDev==TRUE){paste0("logPrevYrSITE.INDEXYearDev+")}else{paste0("logPrevYrSITE.INDEX+")}},
                                    if(length(tempcols)>0){paste0("`",MaxTempVar,"_",WindowModels[focrow,tempcols],"SpatialAv`*`",MaxTempVar,"_",WindowModels[focrow,tempcols],"YearDev`+I(`",MaxTempVar,"_",WindowModels[focrow,tempcols],"SpatialAv`^2)",collapse="",sep="+")},                        
                                    if(length(precipcols)>0){paste0("`",PrecipVar,"_",WindowModels[focrow, precipcols], "SpatialAv`*`",PrecipVar,"_",WindowModels[focrow, precipcols],"YearDev`+I(`",PrecipVar,"_",WindowModels[focrow, precipcols],"SpatialAv`^2)",sep = "+",collapse="")},
                                    "(1|Year)+(1|Site)+(1|GridCell50km)+(1|GridCell50km:Year)+(1|resid)"))
      } 
      
      if(SpaceOnly==TRUE){
        focformula<-as.formula(paste0("SITE.INDEX~",if(AR1==TRUE){if(AR1YearDev==TRUE){paste0("logPrevYrSITE.INDEXYearDev+")}else{paste0("logPrevYrSITE.INDEX+")}},
                                      if(length(tempcols)>0){paste0("`",MaxTempVar,"_",WindowModels[focrow,tempcols],"SpatialAv`+I(`",MaxTempVar,"_",WindowModels[focrow,tempcols],"SpatialAv`^2)",collapse="",sep="+")},                        
                                      if(length(precipcols)>0){paste0("`",PrecipVar,"_",WindowModels[focrow, precipcols], "SpatialAv`+I(`",PrecipVar,"_",WindowModels[focrow, precipcols],"SpatialAv`^2)",sep = "+",collapse="")},
                                      "(1|Year)+(1|Site)+(1|GridCell50km)+(1|GridCell50km:Year)+(1|resid)"))
      }
      
      if(TimeOnly==TRUE){
        focformula<-as.formula(paste0("SITE.INDEX~",if(AR1==TRUE){if(AR1YearDev==TRUE){paste0("logPrevYrSITE.INDEXYearDev+")}else{paste0("logPrevYrSITE.INDEX+")}},
                                      if(length(tempcols)>0){paste0("`",MaxTempVar,"_",WindowModels[focrow,tempcols],"YearDev`+`",MaxTempVar,"_",WindowModels[focrow,tempcols],"YearDev`:`",MaxTempVar,"_",WindowModels[focrow,tempcols],"SpatialAv`",collapse="",sep="+")},                        
                                      if(length(precipcols)>0){paste0("`",PrecipVar,"_",WindowModels[focrow, precipcols],"YearDev`+`",PrecipVar,"_",WindowModels[focrow, precipcols], "YearDev`:`",PrecipVar,"_",WindowModels[focrow, precipcols],"SpatialAv`",sep = "+",collapse="")},
                                      "(1|Year)+(1|Site)+(1|GridCell50km)+(1|GridCell50km:Year)+(1|resid)"))
      }
      

      # if using Site Indices, we want a Poisson model
      focmod<-try(glmer(focformula,family=poisson,data=AbundanceData,na.action=na.omit,start=if(StartVals==TRUE){list(fixef=StartValues)}else{NULL}))
      
      # if there was an error, try another optimiser
      if(class(focmod)=="try-error"){
        
        print("fitting error, trying another optimiser")
        focmod<-try(glmer(focformula,family=poisson,data=AbundanceData,na.action=na.omit,glmerControl(optimizer="bobyqa"),start=if(StartVals==TRUE){list(fixef=StartValues)}else{NULL}))
      
        }else{
        
      # or if there is a convergence warning, try a different optimiser
      if(length(grep("Model failed to converge with max",summary(focmod)$optinfo$conv$lme4$messages))==1){
        print("failed to converge, trying another optimiser")
        focmod<-try(glmer(focformula,family=poisson,data=AbundanceData,na.action=na.omit,glmerControl(optimizer="bobyqa"),start=if(StartVals==TRUE){list(fixef=StartValues)}else{NULL}))
      }
      }
      
    
    # save model outputs in table 
    
    if(class(focmod)=="try-error"){ 
      print(paste0(focrow," NOT RUN"))
      }else{# if the model ran, input estimates (else this will remain as NAs)
        
    WindowModels$AIC[focrow]<-ModifiedAIC(model=focmod,nwindows=(length(tempcols)+length(precipcols)))
    WindowModels$logLik[focrow]<-logLik(focmod)
    
    if(SpeciesName=="Anthocharis_cardamines"){
      WindowModels$R2M[focrow]<-partR2(focmod,data=AbundanceData,R2_type = "marginal")$R2$estimate
      WindowModels$R2C[focrow]<-partR2(focmod,data=AbundanceData,R2_type = "conditional")$R2$estimate
    }
    
    WindowModels$InterceptEst[focrow]<-summary(focmod)$coef["(Intercept)","Estimate"]
    WindowModels$InterceptSE[focrow]<-summary(focmod)$coef["(Intercept)","Std. Error"]
    
    if(AR1==TRUE){
      if(AR1YearDev==TRUE){
        WindowModels$logPrevYrSITE.INDEXYearDevEst[focrow]<-summary(focmod)$coef["logPrevYrSITE.INDEXYearDev","Estimate"]
        WindowModels$logPrevYrSITE.INDEXYearDevSE[focrow]<-summary(focmod)$coef["logPrevYrSITE.INDEXYearDev","Std. Error"]
      }else{
      WindowModels$logPrevYrSITE.INDEXEst[focrow]<-summary(focmod)$coef["logPrevYrSITE.INDEX","Estimate"]
      WindowModels$logPrevYrSITE.INDEXSE[focrow]<-summary(focmod)$coef["logPrevYrSITE.INDEX","Std. Error"]
      }
    }
    
    # estimates for focal window
    if(TimeOnly!=TRUE){
    WindowModels[focrow,paste0(FocWindow,"SpatialAvEst")]<-summary(focmod)$coef[paste0("`",focVar,"_",WindowModels[focrow, paste0("Window",FocWindow)],"SpatialAv`"),"Estimate"]
    WindowModels[focrow,paste0(FocWindow,"SpatialAvSE")]<-summary(focmod)$coef[paste0("`",focVar,"_",WindowModels[focrow, paste0("Window",FocWindow)],"SpatialAv`"),"Std. Error"]
    
    WindowModels[focrow,paste0(FocWindow,"SpatialAvSqEst")]<-summary(focmod)$coef[paste0("I(`",focVar,"_",WindowModels[focrow, paste0("Window",FocWindow)],"SpatialAv`^2)"),"Estimate"]
    WindowModels[focrow,paste0(FocWindow,"SpatialAvSqSE")]<-summary(focmod)$coef[paste0("I(`",focVar,"_",WindowModels[focrow, paste0("Window",FocWindow)],"SpatialAv`^2)"),"Std. Error"]
    }  
    
    if(SpaceOnly!=TRUE){
    WindowModels[focrow,paste0(FocWindow,"YearDevEst")]<-summary(focmod)$coef[paste0("`",focVar,"_",WindowModels[focrow, paste0("Window",FocWindow)],"YearDev`"),"Estimate"]
    WindowModels[focrow,paste0(FocWindow,"YearDevSE")]<-summary(focmod)$coef[paste0("`",focVar,"_",WindowModels[focrow, paste0("Window",FocWindow)],"YearDev`"),"Std. Error"]
    
    if(TimeOnly==TRUE){
      WindowModels[focrow,paste0(FocWindow,"SpatialAvYearDevEst")]<-summary(focmod)$coef[paste0("`",focVar,"_",WindowModels[focrow, paste0("Window",FocWindow)],"YearDev`:`",focVar,"_",WindowModels[focrow, paste0("Window",FocWindow)],"SpatialAv`"),"Estimate"]
      WindowModels[focrow,paste0(FocWindow,"SpatialAvYearDevSE")]<-summary(focmod)$coef[paste0("`",focVar,"_",WindowModels[focrow, paste0("Window",FocWindow)],"YearDev`:`",focVar,"_",WindowModels[focrow, paste0("Window",FocWindow)],"SpatialAv`"),"Std. Error"]
      
    }else{
      WindowModels[focrow,paste0(FocWindow,"SpatialAvYearDevEst")]<-summary(focmod)$coef[paste0("`",focVar,"_",WindowModels[focrow, paste0("Window",FocWindow)],"SpatialAv`:`",focVar,"_",WindowModels[focrow, paste0("Window",FocWindow)],"YearDev`"),"Estimate"]
      WindowModels[focrow,paste0(FocWindow,"SpatialAvYearDevSE")]<-summary(focmod)$coef[paste0("`",focVar,"_",WindowModels[focrow, paste0("Window",FocWindow)],"SpatialAv`:`",focVar,"_",WindowModels[focrow, paste0("Window",FocWindow)],"YearDev`"),"Std. Error"]
    }
    }
    
    # estimates for previous windows 
    if(length(PrevWindows)>0){
      for(prev in PrevWindows){
        
        # is this window temperature or precipitation?
        # this is needed for extracting model estimates 
        if(length(grep("MaxTemp",prev))==1){PrevVar<-MaxTempVar}
        if(length(grep("Precip",prev))==1){PrevVar<-PrecipVar}
        
        if(TimeOnly!=TRUE){
        WindowModels[focrow,paste0(prev,"SpatialAvEst")]<-summary(focmod)$coef[paste0("`",PrevVar,"_",WindowModels[focrow, paste0("Window",prev)],"SpatialAv`"),"Estimate"]
        WindowModels[focrow,paste0(prev,"SpatialAvSE")]<-summary(focmod)$coef[paste0("`",PrevVar,"_",WindowModels[focrow, paste0("Window",prev)],"SpatialAv`"),"Std. Error"]
        
        WindowModels[focrow,paste0(prev,"SpatialAvSqEst")]<-summary(focmod)$coef[paste0("I(`",PrevVar,"_",WindowModels[focrow, paste0("Window",prev)],"SpatialAv`^2)"),"Estimate"]
        WindowModels[focrow,paste0(prev,"SpatialAvSqSE")]<-summary(focmod)$coef[paste0("I(`",PrevVar,"_",WindowModels[focrow, paste0("Window",prev)],"SpatialAv`^2)"),"Std. Error"]
        }
        
        if(SpaceOnly!=TRUE){
        WindowModels[focrow,paste0(prev,"YearDevEst")]<-summary(focmod)$coef[paste0("`",PrevVar,"_",WindowModels[focrow, paste0("Window",prev)],"YearDev`"),"Estimate"]
        WindowModels[focrow,paste0(prev,"YearDevSE")]<-summary(focmod)$coef[paste0("`",PrevVar,"_",WindowModels[focrow, paste0("Window",prev)],"YearDev`"),"Std. Error"]
        
        if(TimeOnly==TRUE){
          WindowModels[focrow,paste0(prev,"SpatialAvYearDevEst")]<-summary(focmod)$coef[paste0("`",PrevVar,"_",WindowModels[focrow, paste0("Window",prev)],"YearDev`:`",PrevVar,"_",WindowModels[focrow, paste0("Window",prev)],"SpatialAv`"),"Estimate"]
          WindowModels[focrow,paste0(prev,"SpatialAvYearDevSE")]<-summary(focmod)$coef[paste0("`",PrevVar,"_",WindowModels[focrow, paste0("Window",prev)],"YearDev`:`",PrevVar,"_",WindowModels[focrow, paste0("Window",prev)],"SpatialAv`"),"Std. Error"]
        }else{
        WindowModels[focrow,paste0(prev,"SpatialAvYearDevEst")]<-summary(focmod)$coef[paste0("`",PrevVar,"_",WindowModels[focrow, paste0("Window",prev)],"SpatialAv`:`",PrevVar,"_",WindowModels[focrow, paste0("Window",prev)],"YearDev`"),"Estimate"]
        WindowModels[focrow,paste0(prev,"SpatialAvYearDevSE")]<-summary(focmod)$coef[paste0("`",PrevVar,"_",WindowModels[focrow, paste0("Window",prev)],"SpatialAv`:`",PrevVar,"_",WindowModels[focrow, paste0("Window",prev)],"YearDev`"),"Std. Error"]
        }
        }
      } 
    }
    
    if(length(grep("Model failed to converge",summary(focmod)$optinfo$conv$lme4$messages))==1){
      WindowModels$ConvergeWarning[focrow]<-"YES"
    }
    if(length(grep("singular",summary(focmod)$optinfo$conv$lme4$messages))==1){
      WindowModels$SingularWarning[focrow]<-"YES"
    }
    if(length(grep("eigenvalue",summary(focmod)$optinfo$conv$lme4$messages))==1){
      WindowModels$EigenvalueWarning[focrow]<-"YES"
    }
    }
      
    # save as it goes along
    write.table(WindowModels[focrow,],
                paste0(modfolder,SpeciesName,"_SlidingWindowRun",FocWindow,modname,if(StartVals==TRUE){paste0("_StartVals")},".txt"),
                append=T,sep="\t",row.names=F,
                col.names=!file.exists(paste0(modfolder,SpeciesName,"_SlidingWindowRun",FocWindow,modname,if(StartVals==TRUE){paste0("_StartVals")},".txt"))) # add col names only in the first save (ie if file doesn't exist)
    }
  }
  #  only when all have run, find best model  
  
  WindowModels<-read.table(paste0(modfolder,SpeciesName,"_SlidingWindowRun",FocWindow,modname,if(StartVals==TRUE){paste0("_StartVals")},".txt"),header = TRUE)

  if(nrow(WindowModels)==nwindows){ # if all models have run (or attempted to run) and saved, find the best
    
    windowDur<-as.character(WindowModels[which(WindowModels$AIC==min(WindowModels$AIC,na.rm=T)),paste0("Window",FocWindow)])
    BestWindows["Window",FocWindow]<-windowDur
    
    # if this is the first window, compare to the null model
    if(LastWindow==FALSE){
      
        # null
        if(AR1==TRUE){
          if(AR1YearDev==TRUE){
            nullformula<-as.formula(paste0("SITE.INDEX~logPrevYrSITE.INDEXYearDev+(1|Year)+(1|Site)+(1|GridCell50km)+(1|GridCell50km:Year)+(1|resid)"))
          }else{
            nullformula<-as.formula(paste0("SITE.INDEX~logPrevYrSITE.INDEX+(1|Year)+(1|Site)+(1|GridCell50km)+(1|GridCell50km:Year)+(1|resid)"))
            }
          }
      if(AR1==FALSE){nullformula<-as.formula(paste0("SITE.INDEX~1+(1|Year)+(1|Site)+(1|GridCell50km)+(1|GridCell50km:Year)+(1|resid)"))}
        nullmod<-glmer(nullformula,family=poisson,data=AbundanceData,na.action=na.omit)
      
      nullmodAIC<-ModifiedAIC(model=nullmod,nwindows=0)
      focAIC<-min(WindowModels$AIC,na.rm=T)
      BestWindows["AIC",FocWindow]<-focAIC
      BestWindows["AICChange",FocWindow]<-nullmodAIC-focAIC
      
      
      # if model is no better than null (AIC at least 2 less)
      # exclude the window
      if((nullmodAIC-focAIC)>2){
        print(paste0("AIC diff = ",nullmodAIC-focAIC))
      }else{
        BestWindows["Window",which(colnames(BestWindows)==FocWindow)]<-"NO WINDOW"
        warning("No Window")
        print("No Window")
        print(paste0("AIC diff = ",nullmodAIC-focAIC))
      }
      
    }else{ # if there is a previous window compare model to that 
      
      focAIC<-min(WindowModels$AIC,na.rm=T)
      BestWindows["AIC",FocWindow]<-focAIC
      
      prevAIC<-as.numeric(BestWindows["AIC",LastWindow])
      BestWindows["AICChange",FocWindow]<-prevAIC-focAIC
      
      # if model is no better than previous model (AIC at least 2 less)
      # exclude the window
      if((prevAIC-focAIC)>2){
        print(paste0("AIC diff = ",prevAIC-focAIC))
      }else{
        BestWindows["Window",which(colnames(BestWindows)==FocWindow)]<-"NO WINDOW"
        warning("No Window")
        print("No Window")
        print(paste0("AIC diff = ",prevAIC-focAIC))
      }
    }
    
    print(BestWindows)
  
    print(paste0("ConvergeWarning n=",length(which(!is.na(WindowModels$ConvergeWarning)))))
    print(paste0("SingularWarning n=",length(which(!is.na(WindowModels$SingularWarning)))))
    print(paste0("EigenvalueWarning n=",length(which(!is.na(WindowModels$EigenvalueWarning)))))
    
    print(paste0(length(which(is.na(WindowModels$logLik))), " not run:", paste(WindowModels[which(is.na(WindowModels$logLik)),paste0("Window",FocWindow)],collapse=", ")))
    
    # save models
    # this will be overwritten with each new climate variable window
    
    save(BestWindows, DayNumbers, AbundanceData, 
          file=paste0(modfolder,SpeciesName,"_SlidingWindowOutputs",modname,if(StartVals==TRUE){paste0("_StartVals")},".rda"))
    
    # this is the output to use in the space v time 
    # this will overwrite previous search if start values are used, but also saving the above means that we keep the start value and non-start value outputs for comparison
    save(BestWindows, DayNumbers, AbundanceData, 
         file=paste0(modfolder,SpeciesName,"_SlidingWindowOutputsToUse",modname,".rda"))
    
  }
}

#########################################
#########################################

# function to find modified AIC with the start and end of windows counted as parameters
ModifiedAIC<-function(model,nwindows){
  logLik<-logLik(model)
  k<-attributes(logLik)$df 
  AIC<-(2*(k+(2*nwindows)))-(2*logLik)
  AIC<-as.numeric(AIC)
  return(AIC)
}
