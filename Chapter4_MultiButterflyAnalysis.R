rm(list=ls())

# load functions 
source("Functions/WindowFunctions.R")
source("Functions/SpaceTimeFunction.R")
# see function scripts for an explanation of the arguments they take

# specify species for use in this chapter
focalSpecies<-c("Anthocharis_cardamines","Aphantopus_hyperantus","Colias_croceus","Thymelicus_sylvestris",
                "Callophrys_rubi", "Argynnis_aglaja", "Pyronia_tithonus")

for(SpeciesName in focalSpecies){
  
  print(SpeciesName)
  
  # load abundance data (UKBMS site indices), see DataPreparation.R
  AbundanceData<-read.csv(paste0("UKBMS/Site indices/ukbmsSiteIndices2021_",SpeciesName,".csv"),row.names = 1)
  # UKBMS site index and site location data (note that the publicly available site location data excludes sensitive sites) is available here https://catalogue.ceh.ac.uk/documents/571a676f-6c32-489b-b7ec-18dcc617a9f1
  
  # specify windows to search
  CountPercentCurr<-99 # window search ends when 99% of individuals observed in the current year
  CountPercentPrev<- -30 # window search began 30 days before first count in the previous year
  WindowDurations<-c(seq(20,150,10)) # different window lengths 
  WindowStartInterval<-10 # start windows at 10 day intervals
  
  #################
  ## Run windows ##
  #################
  
  # Get the climate values for the windows with spatial average calculated across the years 1976-90
  GetWindowsClimates(SpeciesName, AbundanceData, CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev, 
                     WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, SpAvSubset=1976:1990)
 
  ## run the sliding window search for this species
  RunSlidingWindow(SpeciesName, Scaled=TRUE, FocWindow="MaxTemp1",
                   CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev,
                   WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                   AR1=TRUE, AR1YearDev=FALSE, SpAvSubset=1976:1990, StartVals=FALSE)

  
  ###############################
  ## Run space time comparison ##
  ###############################
  
  Hurdle<-TRUE
  nitt<-1000000
  thin=100
  
  ## Run the space v time model, including a predictor for the previous year's abundance, as a yearly deviation 
  # (but this is using the window output from a search where the previous year's abundance predictor was not a yearly deviation (above))

  SpaceTimeModel(SpeciesName, AbundanceMeasure, AR1=TRUE, AR1YearDev=TRUE, Hurdle=Hurdle,
                 CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev, 
                 WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval, 
                 nitt=nitt, thin=thin, SpAvSubset=1976:1990, WindowAR1Diff=FALSE, exclude="Precip")
  
  CompareSpaceTime(SpeciesName, AbundanceMeasure,  AR1=TRUE, AR1YearDev=TRUE, Hurdle=Hurdle,
                   CountPercentCurr=CountPercentCurr, CountPercentPrev=CountPercentPrev, 
                   WindowDurations=WindowDurations, WindowStartInterval=WindowStartInterval,
                   SpAvSubset=1976:1990, nitt=nitt, WindowAR1Diff=FALSE, exclude="Precip")
  
}
