### 2025-09-24_Set58_NewSSMomEval.R ###
#Set58 Evaluation for Manuscript
#Goal: Evaluation using corrected steady-state maternal approximations
#Fixing liquid intake and maternal background

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
myscripts = lapply(filelist, source)

#Name the set, preferably with format Setnumber_SetName. 
#This unique string will be used to pull input data and store outputs

setID = "Set58_NewSSMomEval"

# Read in underlying data 

sourceData = Read.Input.Data(setID)

#Logspline information is not needed for lookup or SS moms

#Run model for all scenarios

Model.All.Scenarios(sourceData)

#Collect serum at selected dates for any evaluations, projections, or distribution fits

Extract.Serum.at.Dates(setID, sourceData)
#Check stability at the selected dates for the projections and logspline fits
#Saves a .csv file and warns you which scenarios and dates are unstable.
#This make take several minutes to run. 

Check.Stability(setID, sourceData, percentileset = c(0.5, 0.75, 0.90, 0.95), error = 0.05)

#Create evaluation data frame for comparisons
Evaluate.Against.Observed.Sera(setID = setID, 
                               sourceData = sourceData)

#Groupings for RMSE calculations

mymetricgroups = list("CenterOnly" = c("Metric_0.5", "Metric_geomean"),
                      "UppermostOnly" = c("Metric_0.95", "Metric_0.9"),
                      "p75Only" = c("Metric_0.75"),
                      "p25Only" = c("Metric_0.25"),
                      "MedandUpper" = c("Metric_0.5", "Metric_geomean", "Metric_0.95", "Metric_0.9"),
                      "MedandAbove" = c("Metric_0.5", "Metric_geomean", "Metric_0.95", "Metric_0.9", "Metric_0.75"))

mychemgroups = list("HighInfo" = c("PFOA", "PFOS", "PFHxS"),
                    "LowInfo" = c("PFNA", "PFDA", "PFHpA"))

mypopgroups = list("Infants" = c("infants"),
                   "Non-Infants" = c("children.hiexp", "children.loexp", "adults.loexp", 
                                     "children", "women", "men", "adults"), 
                   "Adults" = c("adults.loexp", "women", "men", "adults"), 
                   "Children" = c("children.hiexp", "children.loexp", "children"), 
                   "Non-Adults" = c("infants", "children.hiexp", "children.loexp", "children"))

myrmsetab = Calculate.RMSE(setID, sourceData, 
               MetricGroups = mymetricgroups,
               ChemicalGroups = mychemgroups,
               PopulationGroups = mypopgroups)

widermse = myrmsetab %>%
  pivot_wider(names_from = ChemicalGroup, values_from = RMSE)

write.csv(widermse, file = here(setID, paste0("WiderRMSETable_", setID, ".csv")))

