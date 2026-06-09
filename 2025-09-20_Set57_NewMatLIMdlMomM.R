### 2025-09-20_Set57_NewMatLIMdlMomM.R ###
#Set57 - New maternal liquid intake for male receptor with modeled mother
#Goal: Model male population using distribution from Set56 females
#to calculate CGWCs at age 6 months

#Install and/or load packages

packages = c("here", # file location referencing
             "readxl", # read excel files
             "dplyr", # data manipulation
             "tidyr", #data manipulation
             "purrr", # pmap function
             "ggplot2", # plotting
             "matrixStats", # row quantiles
             "lubridate", # extract year from date
             "MASS",  #Fit lognormal distribution to output
             "EnvStats", #Triangular distribution for over age 70 when varying age
             "logspline" ) #Use logspline to estimate density

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE, repos='http://cran.us.r-project.org')
      library(x, character.only = TRUE)
    }
  }
)

#Read in all the script files
filelist = list.files(here("scripts"), recursive = TRUE, full.names = TRUE)
lapply(filelist, source)

#Name the set, preferably with format Setnumber_SetName. 
#This unique string will be used to pull input data and store outputs

setID = "Set57_NewMatLIMdlMomM"


# Read in underlying data 

sourceData = Read.Input.Data(setID)

#Attach the logspline information to the sourceData list
#Necessary to use modeled mothers
maternalID = "Set56_NewMatLISSMomF"
load(here(maternalID, paste0("LogsplineFits_", maternalID, ".RData")))
sourceData[["MatLogsplines"]] = logsplinefits

#Run model for all scenarios

Model.All.Scenarios(sourceData)

#Collect serum at selected dates for any evaluations, projections, or distribution fits

Extract.Serum.at.Dates(setID, sourceData)


#Project the water concentrations or durations based on the serum data
#This also plots the residuals from the projection in the folder ResidualPlots
Project.from.Serum(setID, sourceData)


#Generate data table for dependent variable vs serum quantiles plots

Prep.Projection.Plot.Data(setID, sourceData, addlquants = NA)


#Check stability at the selected dates for the projections and logspline fits
#Saves a .csv file and warns you which scenarios and dates are unstable.
#This make take several minutes to run. 

Check.Stability(setID, sourceData, percentileset = c(0.5, 0.75, 0.9), error = 0.05)


#Prepare data to plot serum over lifetime

Prep.Lifetime.Plot.Data(setID, sourceData, addlquants = NA)



#Collect breast milk concentrations out of results
#object is called allMCs with similar structure to results

fgps = unique(sourceData$Scenario$filegroup)
allMCvars = list()

for (onefgp in fgps) {
  print(paste0("Starting processing for ", onefgp, " at ", format(Sys.time(), "%H:%M:%S")))
  
  #Load the data for the filegroup. We use multiple groups because a laptop may struggle with a very large file
  load(here(setID, paste0("RawData_", setID, "_", onefgp, ".RData")))
  
  partMCvars = lapply(allMCs, function(x) {x$MC_vars})
  allMCvars = c(allMCvars, partMCvars)
}

save(allMCvars, file = here(setID, paste0("allMCvars_", setID, ".RData")) )

#Calculate starting breast milk concentration for each of these.
#Add to the dataframe and save again

allMCvars = lapply(allMCvars, function(x) {mutate(x, MilkStart = MatStart*1000* LactXfer * (1-(BirthLoss*PlacXfer)))})

save(allMCvars, file = here(setID, paste0("allMCvars_", setID, ".RData")) )

milkquants = data.frame()

for (i in 1:length(allMCvars)) {
  milkquants[i,"ScenID"] = names(allMCvars)[i]
  milkquants[i,"Milk50p"] = quantile(allMCvars[[i]]$MilkStart, probs = 0.50)
  milkquants[i,"Milk75p"] = quantile(allMCvars[[i]]$MilkStart, probs = 0.75)
  milkquants[i,"Milk90p"] = quantile(allMCvars[[i]]$MilkStart, probs = 0.90)
}

write.csv(milkquants, file = here(setID, paste0("StartingMilkQuantiles_", setID, ".csv")))

