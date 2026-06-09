ReadMe for the MassDEP PFAS TK Model 
Emma Spady, 2025-11-06 
Goals 
This manual is meant to assist in running the population PFAS TK model developed by the 
O
ice of Research and Standards at the Massachusetts Department of Environmental 
Protection. This document is not intended to stand alone – please refer to the manuscript 
for the purpose and applications of the model.  
Basic Instructions 
Extract the zipped file PFAS_TK_Model to a relevant directory. 
If you are testing the code and wish for it to run more quickly, you may wish to decrease the 
number of iterations in each exposure scenario. These can be changed in the inputdata 
Excel files in the Scenario tab under “niter”. I strongly recommend using the same “niter” 
for all exposure scenarios in a single inputdata file. Otherwise, some of the automatic post
processing will not work correctly.  
Open the project “PFAS_TK_Model_MassDEP_V2.4.Rproj” in RStudio. 
Open the first working script “2025-09-19_Set56_NewMatLISSMomF.R” 
Run the functions starting from the top. The package.check process should facilitate 
installing any packages you don’t already have. The correct setID string is already entered.  
The longest process to run is the function Model.All.Scenarios(sourceData), which may 
take 24 hours if using the original 35,000 iterations. It may run substantially faster on a 
more powerful computer.  
After running all functions in “2025-09-19_Set56_NewMatLISSMomF.R”, check that the 
following files have been produced in the Set56_NewMatLISSMomF directory: - 
Six files called RawData_Set56_NewMatLISSMomF_PFAS_F.RData, with each 
chemical replacing “PFAS” in the file name. - 
SerumatDates_Set56_NewMatLISSMomF.RData 
- - - - - 
LogsplineFits_Set56_NewMatLISSMomF.RData 
ProjPlotQuants_Set56_NewMatLISSMomF.RData 
ProjPlotQuants_Set56_NewMatLISSMomF.csv 
StabilityCheck_Set56_NewMatLISSMomF.csv 
ProjectionOutput_Set56_NewMatLISSMomF.csv 
Restart your RStudio session, then re-open the project “PFAS_TK_Model.Rproj”. Open the                                   
second working script “2025-09-20_Set57_NewMatLIMdlMomM.R”. The correct setID and 
maternalID are already entered. 
Run the functions starting from the top. Depending on your number of iterations, 
Model.All.Scenarios(sourceData) may take many hours.  
After running all functions in “2025-09-20_Set57_NewMatLIMdlMomM.R”, check that the 
following files have been produced in the Set57_NewMatLIMdlMomM directory: - 
Six files called RawData_Set57_NewMatLIMdlMomM_PFAS_M.RData, with each 
chemical replacing “PFAS” in the file name. - - - - - - - 
SerumatDates_Set57_NewMatLIMdlMomM.RData 
ProjPlotQuants_Set57_NewMatLIMdlMomM.RData 
AllDatesQuants_Set57_NewMatLIMdlMomM.RData 
StabilityCheck_Set57_NewMatLIMdlMomM.csv 
ProjectionOutput_Set57_NewMatLIMdlMomM.csv 
ProjPlotQuants_Set57_NewMatLIMdlMomM.csv 
AllDatesQuants_Set57_NewMatLIMdlMomM.csv 
Restart your RStudio session, then re-open the project “PFAS_TK_Model.Rproj”. Open the 
third working script “2025-09-24_Set58_NewSSMomEval.R”. The correct setID is already 
entered. 
Run the functions starting from the top. Depending on your number of iterations, 
Model.All.Scenarios(sourceData) may take many hours.  
After running all functions in “2025-09-24_Set58_NewSSMomEval.R”, check that the 
following files have been produced in the Set58_NewSSMomEval directory: - 
Eleven files called RawData_Set58_NewSSMomEval_StudyLocation.RData, with 
each study location and demographic group replacing “StudyLocation” in the file 
name. - - 
SerumatDates_Set58_NewSSMomEval.RData 
EvalTable_Set58_NewSSMomEval.RData 
- - - 
StabilityCheck_Set57_NewMatLIMdlMomM.csv 
EvalTable_Set58_NewSSMomEval.csv 
RMSEtable_Set58_NewSSMomEval.csv 
General Process 
Fundamentally, the TK model uses information about drinking water exposure scenarios to 
calculate distributions of serum PFAS for each scenario. Serum PFAS distributions are 
calculated at regular intervals over time as the population ages. The calculation is 
performed when parameters change, which is typically monthly below the age of 1 and 
annually beyond this age. The model tracks time by date, in order to reference exposure 
events.  
The model main directory, called PFAS_TK_Model_MassDEP_V2.4, contains an R project, 
the working scripts for each model run, a “scripts” directory for the underlying functions, 
and directories for model runs. Each model run is called a “Set”, with a specified ID, e.g. 
“Set56_NewMatLISSMomF”. This setID must be the name of a folder in the R project 
directory. The setID will generally be appended to model-related files inside that directory. 
Within each set, several “exposure scenarios” can be modeled and analyzed, either 
individually or in relation to one another.  
The inputs and analyses for each set are organized into a spreadsheet called 
“inputdata_setID”. Di
erent tabs manage di
erent parts of the modeling and analysis. The 
role of each tab is described in the “ReadMe” tab, along with the purpose of each column 
or row in the spreadsheet. 
Functions in the “scripts” directory are generally well-commented, but I will give a brief 
overview here. In the order of their use, the functions are: 
Read.Input.Data(): Reads and organizes the Excel file to be compatible with subsequent 
functions. 
Model.All.Scenarios(): Organizes the rows of the inputdata file, each of which represents an 
exposure scenario, and runs them through the TK model process. Each exposure scenario 
has an ID to organize it. The process saves when each exposure scenario is complete. If the 
process is interrupted, re-running Model.All.Scenarios() will continue from exposure 
scenarios without an associated output. However, if an exposure scenario needs to be re
run, delete the corresponding RawData file or list item in the RawData file. Otherwise the 
function will not replace the existing object.  
monte_carlo_vars(): Draws parameters from the defined distributions to prepare for the 
modeling process. It creates a dataframe where each row represents a Monte Carlo 
iteration, and each column a parameter (or set of parameters) that will be used in the TK 
model. It uses a helper function make_rnd(), which reads the parameter entries and 
distribution shape to create random values from the specified distribution.  
calc_Serum_PFAS(): The wrapper function for calculating a single iteration of serum 
changing over time. Produces a list of dataframes. Each list index represents an iteration. 
Each dataframe has rows representing dates and columns representing values such as 
consumed liquid concentration, background serum PFAS, and total serum PFAS at each 
date. This list is saved as the “results” object in the “RawData” files. 
create_frame(): Sets up the dates for measuring serum and changing parameters. Before 
age 1, dates default to once per month aligned with the receptors’ birthday. Dates default 
to the receptors’ birthdays after age 1. Dates also include any time the PFAS concentration 
in water changes during the receptor’s lifetime. The time zoom functionality can generate 
extra dates to consider around a timeframe of interest. The function also generates a 
rounded-down age for eventual matching to liquid intakes. Until age 1, it rounds down to 
the nearest month, then divides by 12. Afterwards, it rounds down to the nearest year.  
pop_Liquid_Intake(): Fills in the liquid intake in mL/(kg*day) at each date, multiplying by the 
liquid intake adjustment factor for the appropriate age group.   
pop_Liquid_Conc(): Fills in water or breast milk PFAS concentrations as appropriate by 
date.   
Water is filled in with a simple matching based on dates. Dates where receptor consumes 
breast milk instead of water are identified using fed_mo. The PFAS concentration in the 
milk is calculated using MilkDecline, LactXfer, and MatStart (appearing in the function as 
momserum).   
pop_Bkgd_Serum(): Fills in background serum PFAS for each date. It also makes the 
background PFAS phase in linearly starting from 0 until age 3.   
calc_Total_Serum(): Uses entries in each date to calculate Total_Serum at the next date. 
The receptor begins with Total_Serum based on placental transfer, with no background. The 
second row needs to be made manually because there isn’t a value for bkgdold before birth 
(so all the placental transfer PFAS is available for elimination). Uses the helper function 
nextserum() to perform the calculation using the integrated rate law for the system. 
Extract.Serum.at.Dates(): Combines all iterations’ serum values into a single dataframe, 
where rows are iterations and columns are scenario IDs with dates. Loads each RawData 
file one by one, to avoid overwhelming small computers. May take several hours if using 
35,000 iterations. 
Downstream analyses are relatively self-explanatory. 
 Fit.Serum.Distributions(): For when model output will be used as the maternal serum 
distribution in another exposure scenario. 
 Evaluate.Against.Observed.Sera(): For comparing modeled serum percentiles against 
observed serum percentiles.  
Project.from.Serum(): For calculating CGWCs and other linear regressions. Connects 
exposure scenario parameters (for example, water concentration) to serum concentration 
at a specified percentile.  
Prep.Projection.Plot.Data(): Saves data in the format required to plot the CGWC projection 
process. 
Prep.Lifetime.Plot.Data(): Saves data in the format required to plot serum percentiles over 
the lifetime of a population. 
Check.Stability(): A process for estimating Monte Carlo error based on principles of 
quantile regression. Helps determine if enough iterations were run.  
Calculate.RMSE(): Calculates RMSE between observed and modeled percentiles. 
