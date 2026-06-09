### 2025-09-19_Set56_NewMatLISSMomF.R ###
#Set56 - New maternal liquid intake for female receptor with steady-state mother
#Goal: Model female population at age 40 and fit distribution 
#to use as the maternal serum of 6 month old males in Set57

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
setID = "Set56_NewMatLISSMomF"

# Read in underlying data 
sourceData = Read.Input.Data(setID)

#Populations with SS mothers do not need the logspline information 
#attached to the sourceData list

#Run model for all scenarios
Model.All.Scenarios(sourceData)

#Collect serum at selected dates for any evaluations, projections, or distribution fits
Extract.Serum.at.Dates(setID, sourceData)


#Fit logspline densities to selected scenarios at a date
Fit.Serum.Distributions(setID, sourceData)

#Project the water concentrations or durations based on the serum data
#This also plots the residuals from the projection in the folder ResidualPlots
Project.from.Serum(setID, sourceData)


#Generate data table for dependent variable vs serum quantiles plots

Prep.Projection.Plot.Data(setID, sourceData, addlquants = NA)


#Check stability at the selected dates for the projections and logspline fits
#Saves a .csv file and warns you which scenarios and dates are unstable.
#This make take several minutes to run. 

Check.Stability(setID, sourceData, 
                percentileset = c(0.5, 0.75, 0.9), error = 0.05)

