
rm(list=ls())

# This chapter uses quantinemo software to run individual-based simulations

# https://www2.unil.ch/popgen/softwares/quantinemo/files/quantinemo.pdf

### install quantinemo R package
# library(devtools)
# install_github("frederic-michaud/RQuantiNemo")

library(RQuantiNemo)
library(viridis)

#################################

### Specify the simulation scenarios to run 

Sims<-data.frame(climatechange="NC",dispersalrate=0.02,selectionintensity=15,carryingcapacity=c("Constant","Varying"),OptSlope=5)

climatechange<-c("NC","CC") # no change and continuous change
dispersal<-c(0.02, 0.2) # dispersal values 
selection<-c(15, 30) # selection intensity 
# carrying capacity - flat or climate determined 

# get desired combinations of the above parameters 
Sims<-rbind(Sims,rbind(expand.grid(climatechange=climatechange,dispersalrate=dispersal,selectionintensity=15,carryingcapacity="Varying",OptSlope=5),
                       expand.grid(climatechange=climatechange,dispersalrate=0.02,selectionintensity=selection,carryingcapacity="Varying",OptSlope=5),
                       expand.grid(climatechange=climatechange,dispersalrate=dispersal,selectionintensity=15,carryingcapacity="Constant",OptSlope=5),
                       expand.grid(climatechange=climatechange,dispersalrate=0.02,selectionintensity=selection,carryingcapacity="Constant",OptSlope=5)))
Sims<-unique(Sims)

# add scenarios to test the required number of burnin generations
Burnin<-Sims
Burnin<-Burnin[-which(Burnin$climatechange=="CC"),] # CC and NC are the same in warmup period (no directional change)
Burnin$Burnin<-TRUE
Sims$Burnin<-FALSE
Sims<-rbind(Burnin,Sims)

# add scenario for weak local adaptation
Sims[nrow(Sims)+1,]<-c("NC",0.5,50,"Varying",5,TRUE)
Sims[nrow(Sims)+1,]<-c("NC",0.5,50,"Constant",5,TRUE)
Sims[nrow(Sims)+1,]<-c("NC",0.5,50,"Varying",5,FALSE)
Sims[nrow(Sims)+1,]<-c("CC",0.5,50,"Varying",5,FALSE)
Sims[nrow(Sims)+1,]<-c("NC",0.5,50,"Constant",5,FALSE)
Sims[nrow(Sims)+1,]<-c("CC",0.5,50,"Constant",5,FALSE)

rownames(Sims)<-1:nrow(Sims)

####################################

# specify annual temperature increase 
AnnTempInc<-0.025 

# loop over simulations with different parameters
for(simnumber in 1:nrow(Sims)){ 
  
  # whether it's burnin or not determines the number of generations and replicates
  if(Sims$Burnin[simnumber]==TRUE){
    nwarmup<-800 # number of warmup generations
    ngens<-nwarmup # total number of generations
    nreps<-4 # number of replicates to run 
  }else{
    nwarmup<-600 # number of warmup generations
    ngens<-nwarmup+50 # total number of generations
    nreps<-100 # number of replicates to run 
  }

  # specify the name of the simulation (for folder and file name when saving)
  SimName<-paste0(Sims$climatechange[simnumber],"_dispersal",unlist(strsplit(as.character(Sims$dispersalrate[simnumber]),"[.]"))[2],
                "_selection", Sims$selectionintensity[simnumber],"_k",Sims$carryingcapacity[simnumber])
  # specific details of this simulation (for file name)
  SimSpecifics<-paste0("_slope",gsub("[.]","-",Sims$OptSlope[simnumber]),
                     if(Sims$Burnin[simnumber]==TRUE){paste0("_BurnIn")})

# Set parameters for this simulation

  # identify this simulation's parameters from table
  focselection<-Sims$selectionintensity[simnumber]
  focdispersal<-Sims$dispersalrate[simnumber]
  
  # get the relationship between temperature and optimum
  focslope<-Sims$OptSlope[simnumber]

  AnnTempSD<-1 # sd of annual temperature variation (ie random fluctuations across years)

for(focrep in 1:nreps){ # for each replicate of this simulation
    
    # if this replicate hasn't been run yet, run it 
    if(!file.exists(paste0("SimulationRuns/",if(Sims$Burnin[simnumber]==TRUE){paste0("BurnIn/")},SimName,"/",SimName,"_Parameters",SimSpecifics,"_",focrep,".rda"))){ 
    
    print(paste0(simnumber, "_",focrep))
    
    # Generate temperature and optima values for this replicate
    
    # temperature data stores
    # temperatures at each site at time 1 - sequence temperatures varying by 10 degrees, centred on zero
    temperature<-(seq(-5,5,length.out=30)) 
    foctemperaturetable<-data.frame(site = 1:30,TemperatureT1 = temperature)
    foctemperaturetable[,paste0("TemperatureT",2:ngens)]<-NA # store for other years 
    
    # optima data stores
    focoptimumtable<-data.frame(site = 1:30,optimumT1 = NA)
    focoptimumtable[,paste0("optimumT",2:ngens)]<-NA 
    
  # if there's no directional climate change in this simulation 
  if(Sims$climatechange[simnumber]=="NC"){ 
    
    # find optimum at time 1
    focoptimumtable$optimumT1 <- (foctemperaturetable$TemperatureT1*focslope)  
    
    for(time in 2:ngens){ # for other time points
      
      # simulate each year's temperature 
      # we still have fluctuations in period with no climate change, so this is temperature at time 1 plus a annual deviate 
      foctemperaturetable[,paste0("TemperatureT",time)]<- foctemperaturetable[,"TemperatureT1"] + rnorm(1,0,AnnTempSD)
      
      # find the corresponding optima
      focoptimumtable[,paste0("optimumT",time)] <- (foctemperaturetable[,paste0("TemperatureT",time)]*focslope)  
      
    }
  }
  
  # if there's climate change in this simulation  
  if(Sims$climatechange[simnumber]=="CC"){ 
    
    # find optimum at time 1
    focoptimumtable$optimumT1 <- (foctemperaturetable$TemperatureT1*focslope)  

    for(time in 2:nwarmup){ # years to settle to equilibrium with no climate change
      
      # simulate each year's temperature 
      # we still have fluctuations in period with no climate change, so this is temperature at time 1 plus a annaul deviate 
      foctemperaturetable[,paste0("TemperatureT",time)]<- foctemperaturetable[,"TemperatureT1"] + rnorm(1,0,AnnTempSD)
      
      # find corresponding optima
      focoptimumtable[,paste0("optimumT",time)] <- (foctemperaturetable[,paste0("TemperatureT",time)]*focslope)  
      
    }
    
    for(time in (nwarmup+1):ngens){ # specify optima under climate change 
      
      # global increase in this year - annual increase plus random yearly fluctuation  
      # this is relative to year 1 so that other year's fluctations aren't included (no autocorrelation)
      globalchange<-AnnTempInc*(time-nwarmup) + rnorm(1,0,AnnTempSD)
      
      # Add global change to the previous year's site temperatures
      foctemperaturetable[,paste0("TemperatureT",time)]<- foctemperaturetable[,"TemperatureT1"] + globalchange 
      
      # find optima
      focoptimumtable[,paste0("optimumT",time)]<- foctemperaturetable[,paste0("TemperatureT",time)]*focslope  
      
    }
    
  }
  
  # get optima into the correct format for quantinemo (see quantinemo manual for details)
  # this is (Generation {{Patch1Value Patch2Value.....Patch30Value}}, Generation {{Patch1Value Patch2Value.....Patch30Value}}....)
  focoptima<-c()
  for(i in 1:ngens){focoptima<-c(focoptima,paste0("{", paste0(focoptimumtable[,paste0("optimumT",i)],collapse=" "),"}"))}
  focoptima<-paste(c(1:ngens),focoptima, collapse = ", ") 
  focoptima<-paste0("(",focoptima,")")
  
  # if this simulation's carrying capacity is constant, specify it as 20,000
  if(Sims$carryingcapacity[simnumber]=="Constant"){
                  foccapacity<-20000
                  CCTable<-foccapacity
  }
  # if this simulation's carrying capacity is temperature-determined
  if(Sims$carryingcapacity[simnumber]=="Varying"){ 
    
                # set up store for carrying capacities
                CCTable<-data.frame(site = foctemperaturetable$site)
                CCTable[,paste0("CCT",1:ngens)]<-NA 
    
                # Use quantiles of a normal distribution to generate carrying capacities across sites
                
                # max carrying capacity is at mid temp of 0 - use this to find multiplier to get from dnorm to abundance 
                multiplier<-20000/dnorm(0,0,3)
                
                # find carrying capacities at each temperature using this
                focCC<-c()
                 for(focT in 1:ngens){
                  CCTable[,paste0("CCT",focT)]<-dnorm(foctemperaturetable[,paste0("TemperatureT",focT)],0,3)*multiplier
                  
                  # put this generation's carrying capacity into format for quantinemo (see quantinemo manual for details)
                  # this is (Generation {Patch1Value Patch2Value.....Patch30Value}, Generation {Patch1Value Patch2Value.....Patch30Value}....)
                  focCC<-c(focCC,paste0("{", paste0( CCTable[,paste0("CCT",focT)],collapse=" "),"}"))
                }
                
                # put this all together into into the format for quantinemo
                focCC<-paste(c(1:ngens),focCC, collapse = ", ") 
                focCC<-paste0("(",focCC,")")
                
                foccapacity<-focCC
            }

  
# specify parameters for the simulation in a list for input to quantinemo

focparameters = list(
  "replicates" = 1, # n replicates of the simulation

  "mating_system" =  3, # random mating with two sexes 
  
  "generations" = ngens, # n generations
  "patch_capacity" = foccapacity, # the capacity of each patch
  "patch_number" = 30, # n patches in the metapopulation 
  "patch_ini_size" = 3000, # start size for each patch 
  
  "mating_nb_offspring_model" = 9, # n offspring at each gen is logistically regulated with a stochastic element
  "growth_rate" = 0.48, # growth rate r of population  
  "dispersal_rate" = focdispersal, # emigration rate
  "dispersal_model" = 2, # 1D stepping stone model 
  "dispersal_border_model" = 2, # absorbing boundaries (individuals can move past border patches)
  "dispersal_rate_model" = 1, # density dependent dispersal rate 
  
  ## Quantitative traits
  "quanti_loci" = 25, # n quantitative loci per individual per trait 
  "quanti_nb_trait" = 1, # n quantitative traits
  "quanti_all" = 2, # max number of alleles per locus 
  "quanti_ini_allele_model" = 0, # initial allele frequency varies across patches
  "quanti_mutation_rate" = 2*10^-4, # mutation rate per locus per generation 
  "quanti_allelic_var" = 0.1, # the variance of the normal distribution from which initial allele frequencies are drawn
  
  "quanti_environmental_model" = 0, # VE specified
  "quanti_heritability" = 27, # VE
  
  # Selection 
  "quanti_selection_model" = 1, # stabilising selection
  "quanti_stab_sel_optima" = focoptima, # the selection optimum for each patch for the trait
  "quanti_stab_sel_intensity" = focselection, # selection intensity, ω. Small value = strong selection, large = weak selection.
  #"patch_stab_sel_intensity_var " # variance of normal distribution by which selection intensity varies at each gen
  "selection_level" = 2, # hard selection 
  
  # specify which statistics to save
  "stat" = "{adlt.nbPops
             adlt.nbInd_p
             fitness
             meanW_p
             varW_p
             quanti
             q.VaW
             q.varA_p
             q.meanG_p
             q.varG_p
             q.meanP_p
             q.varP_p}"  
)

# Run Simulation
  focsim <- new("simulation", parameters = focparameters, 
                sim.name = paste0(SimName,SimSpecifics,"_",focrep), 
                sim.dir = paste0("SimulationRuns/",if(Sims$Burnin[simnumber]==TRUE){paste0("BurnIn/")},SimName,"/"))
  run(focsim)

  save(focparameters,CCTable,foctemperaturetable,focoptimumtable,
       file=paste0("SimulationRuns/",if(Sims$Burnin[simnumber]==TRUE){paste0("BurnIn/")},SimName,"/",SimName,"_Parameters",SimSpecifics,"_",focrep,".rda"))
  
  print(paste0(simnumber,"_",focrep,"_done"))

    }
  }
}
