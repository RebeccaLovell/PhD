rm(list=ls())

# load functions 
source("Functions/WindowFunctions.R")
source("Functions/SpaceTimeFunction.R")
# see function scripts for an explanation of the arguments they take

# specify the orange tip
SpeciesName<-"Anthocharis_cardamines"

# load abundance data (UKBMS site indices), see DataPreparation.R
AbundanceData<-read.csv(paste0("UKBMS/Site indices/ukbmsSiteIndices2021_",SpeciesName,".csv"),row.names = 1)
# UKBMS site index and site location data (note that the publicly available site location data excludes sensitive sites) is available here https://catalogue.ceh.ac.uk/documents/571a676f-6c32-489b-b7ec-18dcc617a9f1

# specify windows to search
CountPercentCurr<-99 # sliding window search ends when 99% of individuals observed in the current year
CountPercentPrev<- -30 # sliding window search began 30 days before first count in the previous year
WindowDurations<-c(seq(20,150,10)) # different window lengths 
WindowStartInterval<-10 # start windows at 10 day intervals

#################
## Run windows ##
#################

# Get the climate values for the windows with spatial average calculated across the years 1976-90
GetWindowsClimates(SpeciesName, AbundanceData, CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev, 
                   WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, SpAvSubset=1976:1990)

### Run the sliding window search with the previous year's abundance predictor as a yearly deviation
### This is the main window search model

# run temperature window search 
RunSlidingWindow(SpeciesName, Scaled=TRUE, FocWindow="MaxTemp1",
                 CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                 WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                 AR1=TRUE, AR1YearDev=TRUE, SpAvSubset=1976:1990, StartVals=FALSE)

# run precipitation window search 
RunSlidingWindow(SpeciesName, Scaled=TRUE, FocWindow="Precip1",LastWindow="MaxTemp1", 
                 CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                 WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                 AR1=TRUE, AR1YearDev=TRUE, SpAvSubset=1976:1990, StartVals=FALSE)

################

### Run the sliding window search with the previous year's abundance as a predictor (not as a yearly deviation)

RunSlidingWindow(SpeciesName, Scaled=TRUE, FocWindow="MaxTemp1",
                 CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                 WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                 AR1=TRUE, AR1YearDev=FALSE, SpAvSubset=1976:1990, StartVals=FALSE)

RunSlidingWindow(SpeciesName, Scaled=TRUE, FocWindow="Precip1",LastWindow="MaxTemp1", 
                 CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                 WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                 AR1=TRUE, AR1YearDev=FALSE, SpAvSubset=1976:1990, StartVals=FALSE)

################

### Run the sliding window search with no predictor for the previous year's abundance 

RunSlidingWindow(SpeciesName, Scaled=TRUE, FocWindow="MaxTemp1",
                 CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                 WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                 AR1=FALSE, AR1YearDev=FALSE, SpAvSubset=1976:1990, StartVals=FALSE)

RunSlidingWindow(SpeciesName, Scaled=TRUE, FocWindow="Precip1",LastWindow="MaxTemp1", 
                 CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                 WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                 AR1=FALSE, AR1YearDev=FALSE, SpAvSubset=1976:1990, StartVals=FALSE)

##############################

### Space only the sliding window search 

RunSlidingWindow(SpeciesName, Scaled=TRUE, FocWindow="MaxTemp1",
                 CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                 WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                 AR1=TRUE, AR1YearDev=TRUE, SpAvSubset=1976:1990, StartVals=FALSE, SpaceOnly=TRUE)

RunSlidingWindow(SpeciesName, Scaled=TRUE, FocWindow="Precip1",LastWindow="MaxTemp1", 
                 CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                 WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                 AR1=TRUE, AR1YearDev=TRUE, SpAvSubset=1976:1990, StartVals=FALSE, SpaceOnly=TRUE)

##############################

### Time only the sliding window search 

RunSlidingWindow(SpeciesName, Scaled=TRUE, FocWindow="MaxTemp1",
                 CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                 WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                 AR1=TRUE, AR1YearDev=TRUE, SpAvSubset=1976:1990, StartVals=FALSE, TimeOnly=TRUE)

RunSlidingWindow(SpeciesName, Scaled=TRUE, FocWindow="Precip1",LastWindow="MaxTemp1", 
                 CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                 WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                 AR1=TRUE, AR1YearDev=TRUE, SpAvSubset=1976:1990, StartVals=FALSE, TimeOnly=TRUE)


###############################
## Run space time comparison ##
###############################

Hurdle<-TRUE # hurdle model
nitt<-5000000 # 5 mil iterations
thin<-100 # thinning interval of 100

### Run the space versus time model with a predictor for the previous year's abundance as a yearly deviation
### This is the main model

# run model
SpaceTimeModel(SpeciesName, AR1=TRUE, AR1YearDev=TRUE, Hurdle=Hurdle,
               CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev, 
               WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
               nitt=nitt, thin=thin, SpAvSubset=1976:1990)

# compare effects in space and time 
CompareSpaceTime(SpeciesName, AR1=TRUE, AR1YearDev=TRUE, Hurdle=Hurdle,
                 CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev, 
                 WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval,
                 SpAvSubset=1976:1990, nitt=nitt)

##############################

### Run the space vs time model with a predictor for the previous year's abundance (not as a yearly deviation)

SpaceTimeModel(SpeciesName, AR1=TRUE, AR1YearDev=FALSE, Hurdle=Hurdle,
               CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev, 
               WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
               nitt=nitt, thin=thin, SpAvSubset=1976:1990)

CompareSpaceTime(SpeciesName,  AR1=TRUE, AR1YearDev=FALSE, Hurdle=Hurdle,
                 CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev, 
                 WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval,
                 SpAvSubset=1976:1990, nitt=nitt)

##############################

### Run the space vs time model with no predictor for the previous year's abundance  

SpaceTimeModel(SpeciesName, AR1=FALSE, AR1YearDev=FALSE, Hurdle=Hurdle,
               CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev, 
               WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
               nitt=nitt, thin=thin, SpAvSubset=1976:1990)

CompareSpaceTime(SpeciesName,  AR1=FALSE, AR1YearDev=FALSE, Hurdle=Hurdle,
                 CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev, 
                 WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval,
                 SpAvSubset=1976:1990, nitt=nitt)

##############################

### Run the space vs time model as a non-hurdle model

Hurdle<-FALSE # not a hurdle model
nitt<-2500000 # 2.5 mil iterations
thin=40 # thinning interval of 40 

SpaceTimeModel(SpeciesName, AR1=TRUE, AR1YearDev=TRUE, Hurdle=Hurdle,
               CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev, 
               WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
               nitt=nitt, thin=thin, SpAvSubset=1976:1990)

CompareSpaceTime(SpeciesName,  AR1=TRUE, AR1YearDev=TRUE, Hurdle=Hurdle,
                 CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev, 
                 WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval,
                 SpAvSubset=1976:1990, nitt=nitt)
