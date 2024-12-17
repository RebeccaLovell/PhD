
# this script contains the functions for running a model comparing spatial and temporal climate effects on occupancy and abundance

library(MCMCglmm)

options(warn=1) # print warnings as they occur

#########################################

# Arguments for the below functions for space-time comparison are as follows:
  # SpeciesName - scientific name in the form "Genus_species"
  # Hurdle - TRUE/FALSE whether to fit a hurdle model
  # Scaled - TRUE/FALSE, should climate variables be scaled and centred (z-score standardised)
  # exclude - exclude a variable? "MaxTemp"or "Precip"
  # AR1 - TRUE/FALSE to include previous year's abundance as a predictor
  # AR1YearDev - TRUE/FALSE - is previous year's abundance predictor using the yearly deviation from the site average abundance?
  # CountPercentCurr - the day when window search will end (in the current year), given as a precentage of counts observed
  # CountPercentPrev - the day the window search will start (in the previous year), given as a precentage of counts observed (-30 means 30 days before first count)
  # WindowDurations - the durations of windows to explore
  # WindowStartInterval- the interval between window starts
  # nitt - numebr of iterations to run
  # thin - thinning interval
  # SpAvSubset - the years spatial average should be calculated over (if not specified, all years are used)
  # WindowAR1Diff - TRUE/FALSE - was the previous year's abundance predictor used for the sliding windows different to that for the space vs time model?
  # AR1Windows - TRUE/FALSE - was the previous year's abundance included as a predictor in the sliding window search?
  # AR1YearDevWindows - if AR1Windows=TRUE, was the previous year's abundance a predictor using the yearly deviation from the site average abundance?
  # REsFiltered - TRUE/FALSE - should the random effects that had zero (or very low) variance in the sliding windows models be excluded from the space versus time model?

#########################################


SpaceTimeModel<-function(SpeciesName, Hurdle=FALSE, Scaled=FALSE, exclude=NA, AR1=FALSE,AR1YearDev=FALSE,  
                         CountPercentCurr, CountPercentPrev, WindowDurations, WindowStartInterval, nitt=NA, thin=NA, SpAvSubset=NA,
                         WindowAR1Diff=FALSE, AR1Windows=NA, AR1YearDevWindows=NA,REsFiltered=FALSE){
  
  Vars<-c("MaxTemp", "Precip")
  if(is.na(exclude)==FALSE){
    Vars<-Vars[-(which(Vars==exclude))]
  }
  
  if(WindowAR1Diff==FALSE){
    AR1Windows<-AR1
    AR1YearDevWindows<-AR1YearDev
  }
  
  modfolder<-paste0("SpaceTimeComparison/",SpeciesName,"/",
                    if(AR1==TRUE){paste0("AR1")},
                    if(Hurdle==TRUE){paste0("Hurdle")},"/")
  
  modname<-paste0(if(Scaled==TRUE){paste0("Scaled")}, # scaled?
                  "_",WindowStartInterval,"Int_",WindowDurations[1],"-",WindowDurations[length(WindowDurations)],"Dur_", # start interval and duration
                  paste0(Vars,collapse=""), 
                  if(AR1==TRUE){paste0("_AR1")}, # AR1?
                  if(AR1YearDev==TRUE){paste0("YearDev")},
                  if(!is.na(SpAvSubset[1])){paste0("_",paste0(SpAvSubset[1],"-",SpAvSubset[length(SpAvSubset)]),"SpAv")}, # SpatialAv subset?
                  if(Hurdle==TRUE){paste0("_hu")},# hurdle?
                  if(WindowAR1Diff==TRUE){paste0("_Windows",if(AR1Windows==TRUE){paste0("_AR1")},if(AR1YearDevWindows==TRUE){paste0("YearDev")})}) 

  # load sliding window outputs (this is loading scaled windows run)
  load(paste0("SlidingWindows/",SpeciesName,"/",CountPercentPrev,"-",CountPercentCurr,"/",
              if(AR1Windows==TRUE){paste0("AR1/")}, # AR1?
           SpeciesName,"_SlidingWindowOutputsToUse_Scaled",
           "_",WindowStartInterval,"Int_",WindowDurations[1],"-",WindowDurations[length(WindowDurations)],"Dur",
           if(AR1Windows==TRUE){paste0("_AR1")}, # AR1?
           if(AR1YearDevWindows==TRUE){paste0("YearDev")},
           if(!is.na(SpAvSubset[1])){paste0("_",paste0(SpAvSubset[1],"-",SpAvSubset[length(SpAvSubset)]),"SpAv")},".rda"))
  
  randomef<-c("Site","Year","GridCell50km","GridCell50km:Year") # full set of REs
  
  if(REsFiltered==TRUE){
  # load REs to exclude (those with var zero in window search)
  if(any(Vars=="Precip")){ # if precip or temp and precip model, load the REs to exclude based on the precip window search (this search includes any identified temp window)
    
    load(file=paste0("SlidingWindows/",SpeciesName,"/",CountPercentPrev,"-",CountPercentCurr,"/",
                     if(AR1Windows==TRUE){paste0("AR1/")}, # AR1?
                     SpeciesName,"_REsExclude_Scaled","_",WindowStartInterval,"Int_",WindowDurations[1],"-",
                     WindowDurations[length(WindowDurations)],"Dur",           
                     if(AR1Windows==TRUE){paste0("_AR1")}, # AR1?
                     if(AR1YearDevWindows==TRUE){paste0("YearDev")},
                     if(!is.na(SpAvSubset[1])){paste0("_",paste0(SpAvSubset[1],"-",SpAvSubset[length(SpAvSubset)]),"SpAv")},"Precip1.rda"))

        }else{  # otherwise, load the temperature only one
      load(file=paste0("SlidingWindows/",SpeciesName,"/",CountPercentPrev,"-",CountPercentCurr,"/",
                       if(AR1Windows==TRUE){paste0("AR1/")}, # AR1?
                       SpeciesName,"_REsExclude_Scaled","_",WindowStartInterval,"Int_",WindowDurations[1],"-",
                       WindowDurations[length(WindowDurations)],"Dur",
                       if(AR1Windows==TRUE){paste0("_AR1")}, # AR1?
                       if(AR1YearDevWindows==TRUE){paste0("YearDev")},
                       if(!is.na(SpAvSubset[1])){paste0("_",paste0(SpAvSubset[1],"-",SpAvSubset[length(SpAvSubset)]),"SpAv")},"MaxTemp1.rda"))
    } 
  
    if(length(SDZero)>0){randomef<-randomef[-match(SDZero,randomef)]} # remove these
  }
  
  if(WindowAR1Diff==TRUE && AR1YearDevWindows==FALSE){ # if AR1 yeardev wasn't calculated for windows search 
    AbundanceData$logPrevYrSITE.INDEXSpatialAv<-tapply(AbundanceData$logPrevYrSITE.INDEX,AbundanceData$Site,mean)[as.character(AbundanceData$Site)]
    AbundanceData$logPrevYrSITE.INDEXYearDev<-AbundanceData$logPrevYrSITE.INDEX-AbundanceData$logPrevYrSITE.INDEXSpatialAv
  }
  
  # find columns for variables
  tempcols<-grep("MaxTemp",colnames(BestWindows))
  precipcols<-grep("Precip",colnames(BestWindows))
  
  if(is.na(exclude)==FALSE){ 
  if(exclude=="MaxTemp"){tempcols<-integer()} # if excluding temp, give tempcols length 0
  if(exclude=="Precip"){precipcols<-integer()} # if excluding precip, give precipcols length 0
  
  if(exclude!="MaxTemp"){if(BestWindows["Window",tempcols]=="NO WINDOW"){tempcols<-integer()}} # if no temp window, give tempcols length 0
  if(exclude!="Precip"){if(BestWindows["Window",precipcols]=="NO WINDOW"){precipcols<-integer()}} # if no precip window, give precipcols length 0
  }else{
    
    if(BestWindows["Window",tempcols]=="NO WINDOW"){tempcols<-integer()} # if no temp window, give tempcols length 0
    if(BestWindows["Window",precipcols]=="NO WINDOW"){precipcols<-integer()} # if no precip window, give precipcols length 0
}
  
 
  Windows<-c()
  
  if(length(tempcols)>0){
      Windows<-c(Windows,paste0("MaxTemp",if(Scaled==TRUE){paste0("Scaled")},"_",BestWindows["Window",tempcols]))
  }
  if(length(precipcols)>0){
    Windows<-c(Windows,paste0("Precip",if(Scaled==TRUE){paste0("Scaled")},"_",BestWindows["Window",precipcols]))
  }
  
  if(Hurdle==FALSE){ # if not a hurdle model
  # parameter expanded priors
    a<-1000
  if(length(randomef)>3){
    pa_prior<-list(R=list(V=diag(1), nu=0.002), 
                   G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
                          G2=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
                          G3=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
                          G4=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a)))
    }else{
  pa_prior<-list(R=list(V=diag(1), nu=0.002), 
                G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
                G2=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
                G3=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a)))
    }
  
  # G is random effects, R is residuals, fixed effects are B
  # default priors for fixed effects
  
  # specify model formula to use, based on windows identified
  predictors<-paste0("`",Windows,"SpatialAv`*`",Windows,"YearDev`+I(`",Windows,"SpatialAv`^2)") # adding `` as var names have numbers in
  
  if(AR1==TRUE){
    if(AR1YearDev==TRUE){
      fixedeffs<-as.formula(paste0("SITE.INDEX~logPrevYrSITE.INDEXYearDev+", paste0(predictors,collapse="+")))
    }else{
    fixedeffs<-as.formula(paste0("SITE.INDEX~logPrevYrSITE.INDEX+", paste0(predictors,collapse="+")))
    }
  }else{
    fixedeffs<-as.formula(paste0("SITE.INDEX~", paste0(predictors,collapse="+")))
  }
    spacetimemod<-MCMCglmm(fixed=fixedeffs, random = as.formula(paste0("~ ",paste0(randomef,collapse = "+"))), 
                            family="poisson", data = AbundanceData, nitt=nitt, thin=thin, prior=pa_prior, pr=TRUE)
    # MCMCglmm automatically includes residual error so don't need to add this   
  
  }
  
  if(Hurdle==TRUE){ # hurdle poisson model
      
      # parameter expanded priors
      
      a<-1000
      if(length(randomef)==4){
      pa_prior<-list(R=list(V=diag(2), nu=0.002, fix=2), # binomial bit is binary, cant est resid have to fix it  
                     G=list(G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a),
                            G2=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a),
                            G3=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a),
                            G4=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a)))
      }else{
        
        if(length(randomef)==3){
        pa_prior<-list(R=list(V=diag(2), nu=0.002, fix=2), # binomial bit is binary, cant est resid have to fix it  
                       G=list(G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a),
                              G2=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a),
                              G3=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a)))
        }else{
          if(length(randomef)==2){
          pa_prior<-list(R=list(V=diag(2), nu=0.002, fix=2), # binomial bit is binary, cant est resid have to fix it  
                         G=list(G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a),
                                G2=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a)))
          }else{
              pa_prior<-list(R=list(V=diag(2), nu=0.002, fix=2), # binomial bit is binary, cant est resid have to fix it  
                             G=list(G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a)))
           }
        }
      }
      
      # G is random effects, R is residuals, fixed effects are B
      # default priors for fixed effects
      
      # specify model formula to use, based on windows identified
      predictors<-paste0("trait:`",Windows,"SpatialAv`*trait:`",Windows,"YearDev`+trait:I(`",Windows,"SpatialAv`^2)") # adding `` as var names have numbers in
      
      if(AR1==TRUE){
        if(AR1YearDev==TRUE){
          fixedeffs<-as.formula(paste0("SITE.INDEX~trait-1 +at.level(trait,1):logPrevYrSITE.INDEXYearDev+", paste0(predictors,collapse="+")))
          #variable:at.level(trait,1) means that the predictor is only in part of the model (in this case the poisson)
          
        }else{
        fixedeffs<-as.formula(paste0("SITE.INDEX~trait-1 +at.level(trait,1):logPrevYrSITE.INDEX+", paste0(predictors,collapse="+")))
        #variable:at.level(trait,1) means that the predictor is only in part of the model (in this case the poisson)
        
        }
      }else{
        fixedeffs<-as.formula(paste0("SITE.INDEX~trait-1 +", paste0(predictors,collapse="+"))) 
      }
      
        spacetimemod<-MCMCglmm(fixed=fixedeffs, random = as.formula(paste0("~ ", paste0(paste0("us(trait):",randomef) ,collapse="+"))), rcov=~idh(trait):units, 
                               family="hupoisson", data = AbundanceData, nitt=nitt, thin=thin, prior=pa_prior, pr=TRUE)
        # us = unstructured - allows covar between binom and pois parts of the model
        
    }
  
  
  #################################
  
  save(spacetimemod, Windows, DayNumbers, 
       file=paste0(modfolder,SpeciesName,"_SpaceTimeModel",modname,"_",nitt,"itt.rda"))

  print(summary(spacetimemod))
}


#########################################
#########################################

CompareSpaceTime<-function(SpeciesName, Hurdle=FALSE, Scaled=FALSE, exclude=NA, AR1=FALSE, AR1YearDev=FALSE,
                           CountPercentCurr, CountPercentPrev, WindowDurations, WindowStartInterval, SpAvSubset=NA, nitt=NA,
                           WindowAR1Diff=FALSE, AR1Windows=NA, AR1YearDevWindows=NA){
  
  Vars<-c("MaxTemp", "Precip")
  if(is.na(exclude)==FALSE){
    Vars<-Vars[-(which(Vars==exclude))]
  }
  
  if(WindowAR1Diff==FALSE){
    AR1Windows<-AR1
    AR1YearDevWindows<-AR1YearDev
  }
  
  # load windows searched 
  load(paste0("SlidingWindows/",SpeciesName,"/", CountPercentPrev,"-",CountPercentCurr,"/",SpeciesName,"_",
              "Windows_",WindowStartInterval,"Int_",
              WindowDurations[1],"-",WindowDurations[length(WindowDurations)],"Dur",
              if(!is.na(SpAvSubset[1])){paste0("_",paste0(SpAvSubset[1],"-",SpAvSubset[length(SpAvSubset)]),"SpAv")},".rda"))
  
  # load sliding window outputs and abundance data (this is loading scaled windows run)
  load(paste0("SlidingWindows/",SpeciesName,"/",CountPercentPrev,"-",CountPercentCurr,"/",
              if(AR1Windows==TRUE){paste0("AR1/")}, # AR1?
              SpeciesName,"_SlidingWindowOutputsToUse_Scaled",
              "_",WindowStartInterval,"Int_",WindowDurations[1],"-",WindowDurations[length(WindowDurations)],"Dur",
              if(AR1Windows==TRUE){paste0("_AR1")}, # AR1?
              if(AR1YearDevWindows==TRUE){paste0("YearDev")},
              if(!is.na(SpAvSubset[1])){paste0("_",paste0(SpAvSubset[1],"-",SpAvSubset[length(SpAvSubset)]),"SpAv")},".rda"))
  
  # load space v time 
  
  modfolder<-paste0("SpaceTimeComparison/",SpeciesName,"/",
                    if(AR1==TRUE){paste0("AR1")},
                    if(Hurdle==TRUE){paste0("Hurdle")},"/")
  
  modname<-paste0(if(Scaled==TRUE){paste0("Scaled")}, # scaled?
                  "_",WindowStartInterval,"Int_",WindowDurations[1],"-",WindowDurations[length(WindowDurations)],"Dur_", # start interval and duration
                  paste0(Vars,collapse=""), 
                  if(AR1==TRUE){paste0("_AR1")}, # AR1?
                  if(AR1YearDev==TRUE){paste0("YearDev")},
                  if(!is.na(SpAvSubset[1])){paste0("_",paste0(SpAvSubset[1],"-",SpAvSubset[length(SpAvSubset)]),"SpAv")}, # SpatialAv subset?
                  if(Hurdle==TRUE){paste0("_hu")},# hurdle?
                  if(WindowAR1Diff==TRUE){paste0("_Windows",if(AR1Windows==TRUE){paste0("_AR1")},if(AR1YearDevWindows==TRUE){paste0("YearDev")})}) 
  
  
  # load space-time comparison model output
  load(paste0(modfolder,SpeciesName,"_SpaceTimeModel",modname,"_",nitt,"itt.rda"))
  
  Optimums<-as.data.frame(matrix(nrow=length(Windows),ncol=2))
  colnames(Optimums)<-c("SpatialOptimum","TemporalOptimum")
  rownames(Optimums)<-Windows
  if(Hurdle==TRUE){ # if hurdle model
    Optimums$SpatialOptimumBinom<-NA  
    Optimums$TemporalOptimumBinom<-NA  
  }
  
  # find the optimum climate in space and in time for each window 
  
  for(i in 1:length(Windows)){ # for each window
    
    focwindow<-Windows[i]
    focrow<-which(rownames(Optimums)==focwindow)
    
    if(Hurdle==FALSE){ # if not a hurdle model
      
      # spatial optimum
      # slope of spatial curve at x is f'(x)=2ax+b 
      # if slope is 0, 0=2ax+b, so x=(-b)/(2a) - this is the peak of the curve, and hence the optimum
      # this is finding the average of the peaks of each iteration, so may not match plots
      Optimums$SpatialOptimum[focrow]<-mean((-(spacetimemod$Sol[,paste0("`",focwindow,"SpatialAv`")]))/(2*(spacetimemod$Sol[,paste0("I(`",focwindow,"SpatialAv`^2)")]))) 
      
      # temporal optimum
      # temporal slope is YeardDev + (YearDev:SpatialAv*x)
      # if slope is 0, 0 = YeardDev + (YearDev:SpatialAv*x), so x=-YearDev/YearDev:SpatialAv
      # this is finding the average of the peaks of each iteration, so may not match plots
      Optimums$TemporalOptimum[focrow]<- mean((-(spacetimemod$Sol[,paste0("`",focwindow,"YearDev`")]))/((spacetimemod$Sol[,paste0("`",focwindow,"SpatialAv`:`",focwindow,"YearDev`")])))
    }
    
    if(Hurdle==TRUE){ # if hurdle model
      
      # find optimum from truncated poisson part of model 
      
      # spatial optimum
      # slope of spatial curve at x is f'(x)=2ax+b 
      # if slope is 0, 0=2ax+b, so x=(-b)/(2a) - this is the peak of the curve, and hence the optimum
      # this is finding the average of the peaks of each iteration, so may not match plots
      Optimums$SpatialOptimum[focrow]<-mean((-(spacetimemod$Sol[,paste0("traitSITE.INDEX:`",focwindow,"SpatialAv`")]))/(2*(spacetimemod$Sol[,paste0("traitSITE.INDEX:I(`",focwindow,"SpatialAv`^2)")])))
      Optimums$SpatialOptimumBinom[focrow]<-mean((-(spacetimemod$Sol[,paste0("traithu_SITE.INDEX:`",focwindow,"SpatialAv`")]))/(2*(spacetimemod$Sol[,paste0("traithu_SITE.INDEX:I(`",focwindow,"SpatialAv`^2)")])))
      
      # temporal optimum
      # temporal slope is YeardDev + (YearDev:SpatialAv*x)
      # if slope is 0, 0 = YeardDev + (YearDev:SpatialAv*x), so x=-YearDev/YearDev:SpatialAv
      # this is finding the average of the peaks of each iteration, so may not match plots
      Optimums$TemporalOptimum[focrow]<- mean((-(spacetimemod$Sol[,paste0("traitSITE.INDEX:`",focwindow,"YearDev`")]))/((spacetimemod$Sol[,paste0("traitSITE.INDEX:`",focwindow,"SpatialAv`:`",focwindow,"YearDev`")])))
      Optimums$TemporalOptimumBinom[focrow]<- mean((-(spacetimemod$Sol[,paste0("traithu_SITE.INDEX:`",focwindow,"YearDev`")]))/((spacetimemod$Sol[,paste0("traithu_SITE.INDEX:`",focwindow,"SpatialAv`:`",focwindow,"YearDev`")])))
    }
  }
  print(Optimums)
  
  # plot model 
  
  if(Hurdle==FALSE){ # if not hurdle model
    
    pdf(file=paste0(modfolder,SpeciesName,"_SpaceTimeModelPlot",modname,"_",nitt,"itt.pdf"),width = 10, height = (5*ceiling(length(Vars)/2)))
    
    par(mfrow=c(ceiling(length(Vars)/2),2))
    
    for(i in 1:length(Windows)){
      
      focwindow<-Windows[i]
      OtherWindow<-Windows[-i] 
      focrow<-which(rownames(Optimums)==focwindow)
      
      focwindowdate<-paste0(ClimWindows$EndDate[which(ClimWindows$window==strsplit(focwindow,"_")[[1]][2])], "-",
                            ClimWindows$StartDate[which(ClimWindows$window==strsplit(focwindow,"_")[[1]][2])]," ",
                            ClimWindows$Year[which(ClimWindows$window==strsplit(focwindow,"_")[[1]][2])], "Year")
      
      min<-round(min(AbundanceData[,paste0(focwindow,"SpatialAv")]),1)
      max<-round(max(AbundanceData[,paste0(focwindow,"SpatialAv")]),1)
      spatialdummy<-seq(min,max,0.1)
      yeardevdummy<-seq(-0.25,0.25,0.1)
      spatialvals<-seq(min,max,0.2)
      
      # values extracted from the model are those when others are held at zero
      # want to plot when other spatialavs are at their means
      # so, find the effect at the mean climate value and add it to the response
      # this is just spatialav as yeardev is just a deviation from this 
      
      if(length(OtherWindow)>0){
        OtherWindowmeanclim<-mean(AbundanceData[,paste0(OtherWindow,"SpatialAv")])
      }
      
      # get spatial response
      spatialresponse<-c()
      
      for(SClim in spatialdummy){
        focspatialresponse<-mean(spacetimemod$Sol[,"(Intercept)"]+(spacetimemod$Sol[,paste0("`",focwindow,"SpatialAv`")]*SClim)+
                                   (spacetimemod$Sol[,paste0("I(`",focwindow,"SpatialAv`^2)")]*SClim^2) + 
                                   if(length(OtherWindow)>0){(spacetimemod$Sol[,paste0("`",OtherWindow,"SpatialAv`")]*OtherWindowmeanclim)+(spacetimemod$Sol[,paste0("I(`",OtherWindow,"SpatialAv`^2)")]*OtherWindowmeanclim^2)}else{0})
        
        spatialresponse<-c(spatialresponse,focspatialresponse)
        
      }
      
      plot(spatialdummy, spatialresponse, type="l", main=paste0(focwindow," (",focwindowdate,")"),
           xlab=paste(strsplit(focwindow,"_")[[1]][1],focwindowdate), ylab="Response")
      
      # find the yeardev response at different spatial values
      
      for(focspatialclim in spatialvals){
        
        # here, included spatialtemp and yeardev effects
        # because you're looking at the deviation from a specific spatialtemp
        # spatialtemp value is held constant while yeardev changes  
        
        yeardevresponse<-c()
        
        for(YClim in yeardevdummy){
          
          focyeardevresponse<-mean(spacetimemod$Sol[,"(Intercept)"]+
                                     (spacetimemod$Sol[,paste0("`",focwindow,"SpatialAv`")]*focspatialclim) +
                                     (spacetimemod$Sol[,paste0("I(`",focwindow,"SpatialAv`^2)")]*focspatialclim^2) +
                                     (spacetimemod$Sol[,paste0("`",focwindow,"YearDev`")]*YClim) +
                                     (spacetimemod$Sol[,paste0("`",focwindow,"SpatialAv`:`",focwindow,"YearDev`")]*focspatialclim*YClim) + 
                                     if(length(OtherWindow)>0){(spacetimemod$Sol[,paste0("`",OtherWindow,"SpatialAv`")]*OtherWindowmeanclim)+(spacetimemod$Sol[,paste0("I(`",OtherWindow,"SpatialAv`^2)")]*OtherWindowmeanclim^2)}else{0})
          
          yeardevresponse<-c(yeardevresponse, focyeardevresponse)
          
        }
        
        # add yeardev line (for this spatialtemp value) to plot  
        # here, yeardevdummy is added to spatialtemp as it's the deviation not the absolute temp value
        lines((focspatialclim+yeardevdummy), yeardevresponse,col="red")  
        
      }
      abline(v=Optimums$SpatialOptimum[focrow], lty=2)
      abline(v=Optimums$TemporalOptimum[focrow],col="red",lty=2)
      mtext(paste0("Optima: Spatial=",round(Optimums$SpatialOptimum[focrow],2),", Temporal=",round(Optimums$TemporalOptimum[focrow],2)),side=3,cex=0.75)
    }
    dev.off()
    
    # Find the posterior of spatial and temporal slopes at different climate values
    # then find the difference between these 
    
    SpaceTimeTable<-data.frame(matrix(nrow=0,ncol=11))
    colnames(SpaceTimeTable)<-c("Window","FocClimVal","SpatialResponseMean","SpatialResponseCredIntLower","SpatialResponseCredIntUpper",
                                "TemporalResponseMean","TemporalResponseCredIntLower", "TemporalResponseCredIntUpper",
                                "DiffMean","DiffCredIntLower","DiffCredIntUpper")
    
    
    pdf(file=paste0(modfolder,SpeciesName,"_SpaceTimeDiffPlot",modname,"_",nitt,"itt.pdf"),width = 15, height = (5*length(Vars)))
    
    par(mfrow=c(length(Vars),3))
    
    for(i in 1:length(Windows)){
      
      focwindow<-Windows[i]
      
      min<-round(min(AbundanceData[,paste0(focwindow,"SpatialAv")]),1)
      max<-round(max(AbundanceData[,paste0(focwindow,"SpatialAv")]),1)
      spatialvals<-seq(min,max,0.1)
      
      # find the rows to fill in
      focrows<-(nrow(SpaceTimeTable)+1):(nrow(SpaceTimeTable)+length(spatialvals))
      SpaceTimeTable[c(focrows),]<-NA
      SpaceTimeTable$Window[focrows]<-focwindow
      SpaceTimeTable$FocClimVal[focrows]<-spatialvals
      
      focwindowdate<-paste0(ClimWindows$EndDate[which(ClimWindows$window==strsplit(focwindow,"_")[[1]][2])], "-",
                            ClimWindows$StartDate[which(ClimWindows$window==strsplit(focwindow,"_")[[1]][2])]," ",
                            ClimWindows$Year[which(ClimWindows$window==strsplit(focwindow,"_")[[1]][2])], " Year")
      
      for(climrow in focrows){ # for each climate value for this window
        focalclim<-SpaceTimeTable$FocClimVal[climrow]
        
        # find the spatial slope posterior at focalclim
        # this is the slope of the tangent that meets the quadratic curve at focalclim
        # to find the slope, find the derivative (using differentiation of the quadratic equation) 
        # if f(x)=ax^2+bx+c, f'(x)=2ax+b - this is the slope of the curve at x
        focspatialslope <- (2*spacetimemod$Sol[,paste0("I(`",focwindow,"SpatialAv`^2)")]*focalclim) + 
          spacetimemod$Sol[,paste0("`",focwindow,"SpatialAv`")]
        SpaceTimeTable$SpatialResponseMean[climrow]<-mean(focspatialslope) # spatial mean
        SpaceTimeTable$SpatialResponseCredIntLower[climrow]<-HPDinterval(focspatialslope)[1,1] # credible interval
        SpaceTimeTable$SpatialResponseCredIntUpper[climrow]<-HPDinterval(focspatialslope)[1,2] # credible interval
        
        # find the temporal posterior at focalclim
        foctemporalslope <- (spacetimemod$Sol[,paste0("`",focwindow,"YearDev`")]) +
          (spacetimemod$Sol[,paste0("`",focwindow,"SpatialAv`:`",focwindow,"YearDev`")]*focalclim) 
        SpaceTimeTable$TemporalResponseMean[climrow]<-mean(foctemporalslope) # temporal mean
        SpaceTimeTable$TemporalResponseCredIntLower[climrow]<-HPDinterval(foctemporalslope)[1,1] # credible interval
        SpaceTimeTable$TemporalResponseCredIntUpper[climrow]<-HPDinterval(foctemporalslope)[1,2] # credible interval
        
        # find the difference between spatial and temporal slopes
        SpatialTemporalDiff<-focspatialslope-foctemporalslope
        SpaceTimeTable$DiffMean[climrow]<-mean(SpatialTemporalDiff) # mean difference
        SpaceTimeTable$DiffCredIntLower[climrow]<-HPDinterval(SpatialTemporalDiff)[1,1] # credible interval of the difference
        SpaceTimeTable$DiffCredIntUpper[climrow]<-HPDinterval(SpatialTemporalDiff)[1,2] # credible interval of the difference
      }
      
      # plot focal clim vs difference mean 
      plot(SpaceTimeTable$FocClimVal[focrows],SpaceTimeTable$DiffMean[focrows], pch=16, main = paste0(focwindow,"Diff"),
           ylim=c(min(SpaceTimeTable$DiffCredIntLower[focrows]),max(SpaceTimeTable$DiffCredIntUpper[focrows])),
           xlab=paste(strsplit(focwindow,"_")[[1]][1],focwindowdate), ylab="Slope Difference (Space-Time)")
      segments(SpaceTimeTable$FocClimVal[focrows],SpaceTimeTable$DiffCredIntLower[focrows],
               SpaceTimeTable$FocClimVal[focrows],SpaceTimeTable$DiffCredIntUpper[focrows], col="red")
      abline(h=0,lty=2)
      dontoverlapzero<-c(SpaceTimeTable$FocClimVal[focrows[which(SpaceTimeTable$DiffCredIntLower[focrows]>0)]],SpaceTimeTable$FocClimVal[focrows[which(SpaceTimeTable$DiffCredIntUpper[focrows]<0)]])
      dontoverlapzero<-dontoverlapzero[order(dontoverlapzero)]
      if(length(dontoverlapzero>0)){
        mtext(paste("diverges from zero at:",paste0(unlist(lapply(split(dontoverlapzero,cumsum(c(0,round(diff(dontoverlapzero),1)!=0.1))),function(x){paste0(range(x)[1],":",range(x)[2])})),collapse=" & ")))
      }
      
      plot(SpaceTimeTable$FocClimVal[focrows],SpaceTimeTable$SpatialResponseMean[focrows], pch=16, main = paste0(focwindow,"SpatialResponse"),
           ylim=c(min(SpaceTimeTable$SpatialResponseCredIntLower[focrows]),max(SpaceTimeTable$SpatialResponseCredIntUpper[focrows])),
           xlab=paste(strsplit(focwindow,"_")[[1]][1],focwindowdate))
      segments(SpaceTimeTable$FocClimVal[focrows],SpaceTimeTable$SpatialResponseCredIntLower[focrows],
               SpaceTimeTable$FocClimVal[focrows],SpaceTimeTable$SpatialResponseCredIntUpper[focrows], col="red")
      abline(h=0,lty=2)
      dontoverlapzero<-c(SpaceTimeTable$FocClimVal[focrows[which(SpaceTimeTable$SpatialResponseCredIntLower[focrows]>0)]],SpaceTimeTable$FocClimVal[focrows[which(SpaceTimeTable$SpatialResponseCredIntUpper[focrows]<0)]])
      dontoverlapzero<-dontoverlapzero[order(dontoverlapzero)]
      if(length(dontoverlapzero>0)){
        mtext(paste("diverges from zero at:",paste0(unlist(lapply(split(dontoverlapzero,cumsum(c(0,round(diff(dontoverlapzero),1)!=0.1))),function(x){paste0(range(x)[1],":",range(x)[2])})),collapse=" & ")))
      }
      
      plot(SpaceTimeTable$FocClimVal[focrows],SpaceTimeTable$TemporalResponseMean[focrows], pch=16, main = paste0(focwindow,"TemporalResponse"),
           ylim=c(min(SpaceTimeTable$TemporalResponseCredIntLower[focrows]),max(SpaceTimeTable$TemporalResponseCredIntUpper[focrows])),
           xlab=paste(strsplit(focwindow,"_")[[1]][1],focwindowdate))
      segments(SpaceTimeTable$FocClimVal[focrows],SpaceTimeTable$TemporalResponseCredIntLower[focrows],
               SpaceTimeTable$FocClimVal[focrows],SpaceTimeTable$TemporalResponseCredIntUpper[focrows], col="red")
      abline(h=0,lty=2)
      dontoverlapzero<-c(SpaceTimeTable$FocClimVal[focrows[which(SpaceTimeTable$TemporalResponseCredIntLower[focrows]>0)]],SpaceTimeTable$FocClimVal[focrows[which(SpaceTimeTable$TemporalResponseCredIntUpper[focrows]<0)]])
      dontoverlapzero<-dontoverlapzero[order(dontoverlapzero)]
      if(length(dontoverlapzero>0)){
        mtext(paste("diverges from zero at:",paste0(unlist(lapply(split(dontoverlapzero,cumsum(c(0,round(diff(dontoverlapzero),1)!=0.1))),function(x){paste0(range(x)[1],":",range(x)[2])})),collapse=" & ")))
      }
      
    }
    
    dev.off()
    
    
    save(SpaceTimeTable, Optimums, Windows, 
         file=paste0(modfolder,SpeciesName,"_SpaceTimeComparison",modname,"_",nitt,"itt.rda"))
    
  }
  
  if(Hurdle==TRUE){ # if hurdle model
    
    pdf(file=paste0(modfolder,SpeciesName,"_SpaceTimeModelPlot",modname,"_",nitt,"itt.pdf"),width = 10, height = (5*length(Vars)))
    
    par(mfrow=c(length(Vars),2))
    
    for(i in 1:length(Windows)){
      
      focwindow<-Windows[i]
      OtherWindow<-Windows[-i] 
      focrow<-which(rownames(Optimums)==focwindow)
      
      focwindowdate<-paste0(ClimWindows$EndDate[which(ClimWindows$window==strsplit(focwindow,"_")[[1]][2])], "-",
                            ClimWindows$StartDate[which(ClimWindows$window==strsplit(focwindow,"_")[[1]][2])]," ",
                            ClimWindows$Year[which(ClimWindows$window==strsplit(focwindow,"_")[[1]][2])], "Year")
      
      min<-round(min(AbundanceData[,paste0(focwindow,"SpatialAv")]),1)
      max<-round(max(AbundanceData[,paste0(focwindow,"SpatialAv")]),1)
      spatialdummy<-seq(min,max,0.1)
      
      yeardevdummy<-seq(-0.25,0.25,0.1)
      spatialvals<-seq(min,max,0.2)
      
      
      # values extracted from the model are those when others are held at zero
      # want to plot when other spatialavs are at their means
      # so, find the effect at the mean climate value and add it to the response
      # this is just spatialav as yeardev is just a deviation from this 
      
      if(length(OtherWindow)>0){
        OtherWindowmeanclim<-mean(AbundanceData[,paste0(OtherWindow,"SpatialAv")])
      }
      
      # get spatial responses
      spatialresponsePois<-c()
      
      # truncated poisson part of model
      for(SClim in spatialdummy){
        
        focspatialresponsePois<- mean((spacetimemod$Sol[,"traitSITE.INDEX"])+(spacetimemod$Sol[,paste0("traitSITE.INDEX:`",focwindow,"SpatialAv`")]*SClim)+
                                        (spacetimemod$Sol[,paste0("traitSITE.INDEX:I(`",focwindow,"SpatialAv`^2)")]*SClim^2)+
                                        if(length(OtherWindow)>0){(spacetimemod$Sol[,paste0("traitSITE.INDEX:`",OtherWindow,"SpatialAv`")]*OtherWindowmeanclim)+(spacetimemod$Sol[,paste0("traitSITE.INDEX:I(`",OtherWindow,"SpatialAv`^2)")]*OtherWindowmeanclim^2)}else{0})
        spatialresponsePois<-c(spatialresponsePois,focspatialresponsePois)
        
      }
      
      plot(spatialdummy, spatialresponsePois, type="l", main=paste0("Pois ",focwindow," (",focwindowdate,")"),
           xlab=paste(strsplit(focwindow,"_")[[1]][1],focwindowdate), ylab="Response")
      
      # find the yeardev response at different spatial values
      
      for(focspatialclim in spatialvals){
        
        # here, included spatialtemp and yeardev effects
        # because you're looking at the deviation from a specific spatialtemp
        # spatialtemp value is held constant while yeardev changes  
        
        yeardevresponsePois<-c()
        
        for(YClim in yeardevdummy){
          
          focyeardevresponsePois<-mean((spacetimemod$Sol[,"traitSITE.INDEX"])+
                                         (spacetimemod$Sol[,paste0("traitSITE.INDEX:`",focwindow,"SpatialAv`")]*focspatialclim) +
                                         (spacetimemod$Sol[,paste0("traitSITE.INDEX:I(`",focwindow,"SpatialAv`^2)")]*focspatialclim^2) +
                                         (spacetimemod$Sol[,paste0("traitSITE.INDEX:`",focwindow,"YearDev`")]*YClim) +
                                         (spacetimemod$Sol[,paste0("traitSITE.INDEX:`",focwindow,"SpatialAv`:`",focwindow,"YearDev`")]*focspatialclim*YClim) + 
                                         if(length(OtherWindow)>0){(spacetimemod$Sol[,paste0("traitSITE.INDEX:`",OtherWindow,"SpatialAv`")]*OtherWindowmeanclim)+
                                             (spacetimemod$Sol[,paste0("traitSITE.INDEX:I(`",OtherWindow,"SpatialAv`^2)")]*OtherWindowmeanclim^2)}else{0})
          
          yeardevresponsePois<-c(yeardevresponsePois, focyeardevresponsePois)
          
        }
        
        # add yeardev line (for this spatialtemp value) to plot  
        # here, yeardevdummy is added to spatialtemp as it's the deviation not the absolute temp value
        lines((focspatialclim+yeardevdummy), yeardevresponsePois,col="red")  
      }
      
      abline(v=Optimums$SpatialOptimum[focrow], lty=2)
      abline(v=Optimums$TemporalOptimum[focrow],col="red",lty=2)
      mtext(paste0("Optima: Spatial=",round(Optimums$SpatialOptimum[focrow],2),", Temporal=",round(Optimums$TemporalOptimum[focrow],2)),side=3,cex=0.75)
      
      # binomial part of model
      spatialresponseBinom<-c()
      
      for(SClim in spatialdummy){
        
        focspatialresponseBinom<-mean((spacetimemod$Sol[,"traithu_SITE.INDEX"])+(spacetimemod$Sol[,paste0("traithu_SITE.INDEX:`",focwindow,"SpatialAv`")]*SClim)+
                                        (spacetimemod$Sol[,paste0("traithu_SITE.INDEX:I(`",focwindow,"SpatialAv`^2)")]*SClim^2) + 
                                        if(length(OtherWindow)>0){(spacetimemod$Sol[,paste0("traithu_SITE.INDEX:`",OtherWindow,"SpatialAv`")]*OtherWindowmeanclim)+
                                            (spacetimemod$Sol[,paste0("traithu_SITE.INDEX:I(`",OtherWindow,"SpatialAv`^2)")]*OtherWindowmeanclim^2)}else{0})
        
        spatialresponseBinom<-c(spatialresponseBinom,focspatialresponseBinom)
        
      }
      
      plot(spatialdummy, spatialresponseBinom, type="l", main=paste0("Binom ",focwindow," (",focwindowdate,")"),
           xlab=paste(strsplit(focwindow,"_")[[1]][1],focwindowdate), ylab="Response")
      
      # find the yeardev response at different spatial values
      
      for(focspatialclim in spatialvals){
        
        # here, included spatialtemp and yeardev effects
        # because you're looking at the deviation from a specific spatialtemp
        # spatialtemp value is held constant while yeardev changes  
        
        yeardevresponseBinom<-c()
        
        for(YClim in yeardevdummy){
          
          focyeardevresponseBinom<-mean((spacetimemod$Sol[,"traithu_SITE.INDEX"])+
                                          (spacetimemod$Sol[,paste0("traithu_SITE.INDEX:`",focwindow,"SpatialAv`")]*focspatialclim) +
                                          (spacetimemod$Sol[,paste0("traithu_SITE.INDEX:I(`",focwindow,"SpatialAv`^2)")]*focspatialclim^2) +
                                          (spacetimemod$Sol[,paste0("traithu_SITE.INDEX:`",focwindow,"YearDev`")]*YClim) +
                                          if(length(OtherWindow)>0){(spacetimemod$Sol[,paste0("traithu_SITE.INDEX:`",focwindow,"SpatialAv`:`",focwindow,"YearDev`")]*focspatialclim*YClim) + 
                                              (spacetimemod$Sol[,paste0("traithu_SITE.INDEX:`",OtherWindow,"SpatialAv`")]*OtherWindowmeanclim)+(spacetimemod$Sol[,paste0("traithu_SITE.INDEX:I(`",OtherWindow,"SpatialAv`^2)")]*OtherWindowmeanclim^2)}else{0})
          
          yeardevresponseBinom<-c(yeardevresponseBinom, focyeardevresponseBinom)
        }
        
        # add yeardev line (for this spatialtemp value) to plot  
        # here, yeardevdummy is added to spatialtemp as it's the deviation not the absolute temp value
        lines((focspatialclim+yeardevdummy), yeardevresponseBinom,col="red")  
      }
      abline(v=Optimums$SpatialOptimumBinom[focrow], lty=2)
      abline(v=Optimums$TemporalOptimumBinom[focrow],col="red",lty=2)
      mtext(paste0("Peak: Spatial=",round(Optimums$SpatialOptimumBinom[focrow],2),", Temporal=",round(Optimums$TemporalOptimumBinom[focrow],2)),side=3,cex=0.75)
      
    }
    dev.off()
    
    # Find the posterior of spatial and temporal slopes at different climate values
    # then find the difference between these 
    
    # first get climate values to explore for each window 
    # these are the whole numbers within the range of SpatialAv
    
    SpaceTimeTable<-data.frame(matrix(nrow=0,ncol=20))
    colnames(SpaceTimeTable)<-c("Window","FocClimVal","SpatialResponseMeanPois","SpatialResponseCredIntLowerPois","SpatialResponseCredIntUpperPois",
                                "TemporalResponseMeanPois","TemporalResponseCredIntLowerPois", "TemporalResponseCredIntUpperPois",
                                "DiffMeanPois","DiffCredIntLowerPois","DiffCredIntUpperPois",
                                "SpatialResponseMeanBinom","SpatialResponseCredIntLowerBinom","SpatialResponseCredIntUpperBinom",
                                "TemporalResponseMeanBinom","TemporalResponseCredIntLowerBinom", "TemporalResponseCredIntUpperBinom",
                                "DiffMeanBinom","DiffCredIntLowerBinom","DiffCredIntUpperBinom")
    
    pdf(file=paste0(modfolder,SpeciesName,"_SpaceTimeDiffPlot",modname,"_",nitt,"itt.pdf"),width = 15, height = (10*(length(Vars))))
    
    par(mfrow=c((2*length(Vars)),3))
    
    for(i in 1:length(Windows)){
      
      focwindow<-Windows[i]
      
      min<-round(min(AbundanceData[,paste0(focwindow,"SpatialAv")]),1)
      max<-round(max(AbundanceData[,paste0(focwindow,"SpatialAv")]),1)
      spatialvals<-seq(min,max,0.1)
      
      # find the rows to fill in
      focrows<-(nrow(SpaceTimeTable)+1):(nrow(SpaceTimeTable)+length(spatialvals))
      SpaceTimeTable[c(focrows),]<-NA
      SpaceTimeTable$Window[focrows]<-focwindow
      SpaceTimeTable$FocClimVal[focrows]<-spatialvals
      
      focwindowdate<-paste0(ClimWindows$EndDate[which(ClimWindows$window==strsplit(focwindow,"_")[[1]][2])], "-",
                            ClimWindows$StartDate[which(ClimWindows$window==strsplit(focwindow,"_")[[1]][2])]," ",
                            ClimWindows$Year[which(ClimWindows$window==strsplit(focwindow,"_")[[1]][2])], " Year")
      
      for(climrow in focrows){ # for each climate value for this window
        focalclim<-SpaceTimeTable$FocClimVal[climrow]
        
        # truncated poisson part 
        
        # find the spatial slope posterior at focalclim
        # this is the slope of the tangent that meets the quadratic curve at focalclim
        # to find the slope, find the derivative (using differentiation of the quadratic equation) 
        # if f(x)=ax^2+bx+c, f'(x)=2ax+b - this is the slope of the curve at x
        focspatialslopePois <- (2*spacetimemod$Sol[,paste0("traitSITE.INDEX:I(`",focwindow,"SpatialAv`^2)")]*focalclim) + 
          spacetimemod$Sol[,paste0("traitSITE.INDEX:`",focwindow,"SpatialAv`")]
        SpaceTimeTable$SpatialResponseMeanPois[climrow]<-mean(focspatialslopePois) # spatial mean
        SpaceTimeTable$SpatialResponseCredIntLowerPois[climrow]<-HPDinterval(focspatialslopePois)[1,1] # credible interval
        SpaceTimeTable$SpatialResponseCredIntUpperPois[climrow]<-HPDinterval(focspatialslopePois)[1,2] # credible interval
        
        # find the temporal posterior at focalclim
        foctemporalslopePois <- (spacetimemod$Sol[,paste0("traitSITE.INDEX:`",focwindow,"YearDev`")]) +
          (spacetimemod$Sol[,paste0("traitSITE.INDEX:`",focwindow,"SpatialAv`:`",focwindow,"YearDev`")]*focalclim) 
        SpaceTimeTable$TemporalResponseMeanPois[climrow]<-mean(foctemporalslopePois) # temporal mean
        SpaceTimeTable$TemporalResponseCredIntLowerPois[climrow]<-HPDinterval(foctemporalslopePois)[1,1] # credible interval
        SpaceTimeTable$TemporalResponseCredIntUpperPois[climrow]<-HPDinterval(foctemporalslopePois)[1,2] # credible interval
        
        # find the difference between spatial and temporal slopes
        SpatialTemporalDiffPois<-focspatialslopePois-foctemporalslopePois
        SpaceTimeTable$DiffMeanPois[climrow]<-mean(SpatialTemporalDiffPois) # mean difference
        SpaceTimeTable$DiffCredIntLowerPois[climrow]<-HPDinterval(SpatialTemporalDiffPois)[1,1] # credible interval of the difference
        SpaceTimeTable$DiffCredIntUpperPois[climrow]<-HPDinterval(SpatialTemporalDiffPois)[1,2] # credible interval of the difference
        
        # binomial part 
        
        # find the spatial slope posterior at focalclim
        # this is the slope of the tangent that meets the quadratic curve at focalclim
        # to find the slope, find the derivative (using differentiation of the quadratic equation) 
        # if f(x)=ax^2+bx+c, f'(x)=2ax+b - this is the slope of the curve at x
        focspatialslopeBinom <- (2*spacetimemod$Sol[,paste0("traithu_SITE.INDEX:I(`",focwindow,"SpatialAv`^2)")]*focalclim) + 
          spacetimemod$Sol[,paste0("traithu_SITE.INDEX:`",focwindow,"SpatialAv`")]
        SpaceTimeTable$SpatialResponseMeanBinom[climrow]<-mean(focspatialslopeBinom) # spatial mean
        SpaceTimeTable$SpatialResponseCredIntLowerBinom[climrow]<-HPDinterval(focspatialslopeBinom)[1,1] # credible interval
        SpaceTimeTable$SpatialResponseCredIntUpperBinom[climrow]<-HPDinterval(focspatialslopeBinom)[1,2] # credible interval
        
        # find the temporal posterior at focalclim
        foctemporalslopeBinom <- (spacetimemod$Sol[,paste0("traithu_SITE.INDEX:`",focwindow,"YearDev`")]) +
          (spacetimemod$Sol[,paste0("traithu_SITE.INDEX:`",focwindow,"SpatialAv`:`",focwindow,"YearDev`")]*focalclim) 
        SpaceTimeTable$TemporalResponseMeanBinom[climrow]<-mean(foctemporalslopeBinom) # temporal mean
        SpaceTimeTable$TemporalResponseCredIntLowerBinom[climrow]<-HPDinterval(foctemporalslopeBinom)[1,1] # credible interval
        SpaceTimeTable$TemporalResponseCredIntUpperBinom[climrow]<-HPDinterval(foctemporalslopeBinom)[1,2] # credible interval
        
        # find the difference between spatial and temporal slopes
        SpatialTemporalDiffBinom<-focspatialslopeBinom-foctemporalslopeBinom
        SpaceTimeTable$DiffMeanBinom[climrow]<-mean(SpatialTemporalDiffBinom) # mean difference
        SpaceTimeTable$DiffCredIntLowerBinom[climrow]<-HPDinterval(SpatialTemporalDiffBinom)[1,1] # credible interval of the difference
        SpaceTimeTable$DiffCredIntUpperBinom[climrow]<-HPDinterval(SpatialTemporalDiffBinom)[1,2] # credible interval of the difference
      }
      
      # plot focal clim vs difference/space/time mean and CI for poisson
      plot(SpaceTimeTable$FocClimVal[focrows],SpaceTimeTable$DiffMeanPois[focrows], pch=16, main = paste("Pois",focwindow),
           ylim=c(min(SpaceTimeTable$DiffCredIntLowerPois[focrows]),max(SpaceTimeTable$DiffCredIntUpperPois[focrows])),
           xlab=paste(strsplit(focwindow,"_")[[1]][1],focwindowdate), ylab="Slope Difference (Space-Time)")
      segments(SpaceTimeTable$FocClimVal[focrows],SpaceTimeTable$DiffCredIntLowerPois[focrows],
               SpaceTimeTable$FocClimVal[focrows],SpaceTimeTable$DiffCredIntUpperPois[focrows], col="red")
      abline(h=0,lty=2)
      dontoverlapzero<-c(SpaceTimeTable$FocClimVal[focrows[which(SpaceTimeTable$DiffCredIntLowerPois[focrows]>0)]],SpaceTimeTable$FocClimVal[focrows[which(SpaceTimeTable$DiffCredIntUpperPois[focrows]<0)]])
      dontoverlapzero<-dontoverlapzero[order(dontoverlapzero)]
      if(length(dontoverlapzero>0)){
        mtext(paste("diverges from zero at:",paste0(unlist(lapply(split(dontoverlapzero,cumsum(c(0,round(diff(dontoverlapzero),1)!=0.1))),function(x){paste0(range(x)[1],":",range(x)[2])})),collapse=" & ")))
      }
      
      plot(SpaceTimeTable$FocClimVal[focrows],SpaceTimeTable$SpatialResponseMeanPois[focrows], pch=16, main = paste0("Pois",focwindow,"SpatialResponse"),
           ylim=c(min(SpaceTimeTable$SpatialResponseCredIntLowerPois[focrows]),max(SpaceTimeTable$SpatialResponseCredIntUpperPois[focrows])),
           xlab=paste(strsplit(focwindow,"_")[[1]][1],focwindowdate))
      segments(SpaceTimeTable$FocClimVal[focrows],SpaceTimeTable$SpatialResponseCredIntLowerPois[focrows],
               SpaceTimeTable$FocClimVal[focrows],SpaceTimeTable$SpatialResponseCredIntUpperPois[focrows], col="red")
      abline(h=0,lty=2)
      dontoverlapzero<-c(SpaceTimeTable$FocClimVal[focrows[which(SpaceTimeTable$SpatialResponseCredIntLowerPois[focrows]>0)]],SpaceTimeTable$FocClimVal[focrows[which(SpaceTimeTable$SpatialResponseCredIntUpperPois[focrows]<0)]])
      dontoverlapzero<-dontoverlapzero[order(dontoverlapzero)]
      if(length(dontoverlapzero>0)){
        mtext(paste("diverges from zero at:",paste0(unlist(lapply(split(dontoverlapzero,cumsum(c(0,round(diff(dontoverlapzero),1)!=0.1))),function(x){paste0(range(x)[1],":",range(x)[2])})),collapse=" & ")))
      }
      
      plot(SpaceTimeTable$FocClimVal[focrows],SpaceTimeTable$TemporalResponseMeanPois[focrows], pch=16, main = paste0("Pois",focwindow,"TemporalResponse"),
           ylim=c(min(SpaceTimeTable$TemporalResponseCredIntLowerPois[focrows]),max(SpaceTimeTable$TemporalResponseCredIntUpperPois[focrows])),
           xlab=paste(strsplit(focwindow,"_")[[1]][1],focwindowdate))
      segments(SpaceTimeTable$FocClimVal[focrows],SpaceTimeTable$TemporalResponseCredIntLowerPois[focrows],
               SpaceTimeTable$FocClimVal[focrows],SpaceTimeTable$TemporalResponseCredIntUpperPois[focrows], col="red")
      abline(h=0,lty=2)
      dontoverlapzero<-c(SpaceTimeTable$FocClimVal[focrows[which(SpaceTimeTable$TemporalResponseCredIntLowerPois[focrows]>0)]],SpaceTimeTable$FocClimVal[focrows[which(SpaceTimeTable$TemporalResponseCredIntUpperPois[focrows]<0)]])
      dontoverlapzero<-dontoverlapzero[order(dontoverlapzero)]
      if(length(dontoverlapzero>0)){
        mtext(paste("diverges from zero at:",paste0(unlist(lapply(split(dontoverlapzero,cumsum(c(0,round(diff(dontoverlapzero),1)!=0.1))),function(x){paste0(range(x)[1],":",range(x)[2])})),collapse=" & ")))
      }
      
      # plot focal clim vs difference/space/time mean and CI for binomial
      plot(SpaceTimeTable$FocClimVal[focrows],SpaceTimeTable$DiffMeanBinom[focrows], pch=16, main = paste("Binom",focwindow),
           ylim=c(min(SpaceTimeTable$DiffCredIntLowerBinom[focrows]),max(SpaceTimeTable$DiffCredIntUpperBinom[focrows])),
           xlab=paste(strsplit(focwindow,"_")[[1]][1],focwindowdate), ylab="Slope Difference (Space-Time)")
      segments(SpaceTimeTable$FocClimVal[focrows],SpaceTimeTable$DiffCredIntLowerBinom[focrows],
               SpaceTimeTable$FocClimVal[focrows],SpaceTimeTable$DiffCredIntUpperBinom[focrows], col="red")
      abline(h=0,lty=2)
      dontoverlapzero<-c(SpaceTimeTable$FocClimVal[focrows[which(SpaceTimeTable$DiffCredIntLowerBinom[focrows]>0)]],SpaceTimeTable$FocClimVal[focrows[which(SpaceTimeTable$DiffCredIntUpperBinom[focrows]<0)]])
      dontoverlapzero<-dontoverlapzero[order(dontoverlapzero)]
      if(length(dontoverlapzero>0)){
        mtext(paste("diverges from zero at:",paste0(unlist(lapply(split(dontoverlapzero,cumsum(c(0,round(diff(dontoverlapzero),1)!=0.1))),function(x){paste0(range(x)[1],":",range(x)[2])})),collapse=" & ")))
      }
      
      plot(SpaceTimeTable$FocClimVal[focrows],SpaceTimeTable$SpatialResponseMeanBinom[focrows], pch=16, main = paste0("Binom",focwindow,"SpatialResponse"),
           ylim=c(min(SpaceTimeTable$SpatialResponseCredIntLowerBinom[focrows]),max(SpaceTimeTable$SpatialResponseCredIntUpperBinom[focrows])),
           xlab=paste(strsplit(focwindow,"_")[[1]][1],focwindowdate))
      segments(SpaceTimeTable$FocClimVal[focrows],SpaceTimeTable$SpatialResponseCredIntLowerBinom[focrows],
               SpaceTimeTable$FocClimVal[focrows],SpaceTimeTable$SpatialResponseCredIntUpperBinom[focrows], col="red")
      abline(h=0,lty=2)
      dontoverlapzero<-c(SpaceTimeTable$FocClimVal[focrows[which(SpaceTimeTable$SpatialResponseCredIntLowerBinom[focrows]>0)]],SpaceTimeTable$FocClimVal[focrows[which(SpaceTimeTable$SpatialResponseCredIntUpperBinom[focrows]<0)]])
      dontoverlapzero<-dontoverlapzero[order(dontoverlapzero)]
      if(length(dontoverlapzero>0)){
        mtext(paste("diverges from zero at:",paste0(unlist(lapply(split(dontoverlapzero,cumsum(c(0,round(diff(dontoverlapzero),1)!=0.1))),function(x){paste0(range(x)[1],":",range(x)[2])})),collapse=" & ")))
      }
      
      plot(SpaceTimeTable$FocClimVal[focrows],SpaceTimeTable$TemporalResponseMeanBinom[focrows], pch=16, main = paste0("Binom",focwindow,"TemporalResponse"),
           ylim=c(min(SpaceTimeTable$TemporalResponseCredIntLowerBinom[focrows]),max(SpaceTimeTable$TemporalResponseCredIntUpperBinom[focrows])),
           xlab=paste(strsplit(focwindow,"_")[[1]][1],focwindowdate))
      segments(SpaceTimeTable$FocClimVal[focrows],SpaceTimeTable$TemporalResponseCredIntLowerBinom[focrows],
               SpaceTimeTable$FocClimVal[focrows],SpaceTimeTable$TemporalResponseCredIntUpperBinom[focrows], col="red")
      abline(h=0,lty=2)
      dontoverlapzero<-c(SpaceTimeTable$FocClimVal[focrows[which(SpaceTimeTable$TemporalResponseCredIntLowerBinom[focrows]>0)]],SpaceTimeTable$FocClimVal[focrows[which(SpaceTimeTable$TemporalResponseCredIntUpperBinom[focrows]<0)]])
      dontoverlapzero<-dontoverlapzero[order(dontoverlapzero)]
      if(length(dontoverlapzero>0)){
        if(length(dontoverlapzero>0)){
          mtext(paste("diverges from zero at:",paste0(unlist(lapply(split(dontoverlapzero,cumsum(c(0,round(diff(dontoverlapzero),1)!=0.1))),function(x){paste0(range(x)[1],":",range(x)[2])})),collapse=" & ")))
        }
      }
      
    }
    
    dev.off()
    
    save(SpaceTimeTable, Optimums, Windows,  
         file=paste0(modfolder,SpeciesName,"_SpaceTimeComparison",modname,"_",nitt,"itt.rda"))
  }
  
}    

