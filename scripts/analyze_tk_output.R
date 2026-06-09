### analyze_tk_output.R ###
#Handles multi-scenario inputs from spreadsheet, as well as the outputs

## Function Calculate.Stable.Iterations() ##
#Estimates how many iterations are needed to be within a particular error at a 
#particular confidence interval for a specified percentile

#ptile = percentile of distribution for evaluation
#error = acceptable error in serum; if absolute, error in ug/L. If percentage, error as number from 0-1. 
#conf = confidence interval for model output, between 0 and 1. I recommend 0.999
#serumvals = vector of serum in ug/L for all iterations

Calculate.Stable.Iterations = function(serumvals, ptile, error, errortype = "absolute", conf) {
  zscore = -qnorm((1-conf)/2, mean = 0, sd = 1)
  ptileval = quantile(serumvals, probs = ptile, na.rm = TRUE)
  
  if (errortype == "percentage") {
    abserror = error*ptileval
  } else if (errortype == "absolute") {
    abserror = error
  }
  sd = abserror/zscore
  
  #Estimate density at desired percentile via logspline 
  #logsplinefit = logspline(serumvals, lbound = 0)
  logsplinefit = Get.Density.Function(serumvals=serumvals)
  estdens = dlogspline(q = ptileval, fit = logsplinefit)
  if (estdens == 0){
    goalniter = 0
  } else {
    goalniter = (ptile*(1-ptile))/((sd*estdens)^2)
  }
  
  return(ceiling(goalniter))
}


## Function Calculate.Percent.Stability() ##
#Estimates the stability achieved at a percentile based on the number of iterations run

Calculate.Percent.Stability = function(serumvals, ptile, niter, conf) {
  
  ptileval = quantile(serumvals, probs = ptile, na.rm = TRUE)
  
  #Estimate density at desired percentile via logspline 
  logsplinefit = Get.Density.Function(serumvals=serumvals)
  estdens = dlogspline(q = ptileval, fit = logsplinefit)
  
  if (estdens == 0) {
    sd = 0
  } else {
    sd = sqrt((ptile*(1-ptile))/(niter*(estdens^2)))
  }
  
  zscore = -qnorm((1-conf)/2, mean = 0, sd = 1)
  
  abserror = sd*zscore
  
  percerror = 100*abserror/ptileval
  
  return(percerror)
}

## Function Get.Density.Function() ##
#Work around errors when using logspline density estimation
#This is the first time ESpady has used try()

Get.Density.Function = function(serumvals, limit.tail = 0.99) {
  dens1 = try(logspline(serumvals, lbound = 0))
  if (class(dens1) == "logspline") {
    print("Density estimate used standard logspline with lbound=0")
    #This is what should usually work
    return(dens1)
  } else {
    dens2 = try(logspline(serumvals))
    if (class(dens2) == "logspline") {
      print("Density estimate required removing lower bound")
      #This was needed when I put together two runs to get 50k iterations in 
      #a CGWC calculation one time and I don't know exactly why it worked
      return(dens2)
    } else {
      truncvals = serumvals
      upperlimit = quantile(serumvals, probs = limit.tail)
      truncvals[which(serumvals>upperlimit)] = upperlimit
      dens3 = try(logspline(truncvals, lbound = 0, ubound = upperlimit))
      if (class(dens3) == "logspline") {
        print("Density estimate required truncation")
        #This is needed if the model output has an extremely long skinny tail
        #which happens with wide distributions of drinking water concentration
        #such as in MaineUS cohort PFOS and PFDA
        return(dens3)
      } else {
        print("No density estimate possible. Returning zero density.")
        print("Stability estimate is not available for this exposure scenario and date.")
        #Aim is for 0s to appear in table but for the stability checker to move to next
        #exposure scenario and date for stability evaluation. 
        dens4 = logspline(runif(length(serumvals),min = -100, max = -5))
        return(dens4)
      }
      
    }
  }
}


## Function Read.Input.Data() ##
#Input manager for Excel spreadsheet
#Converts input spreadsheet to list object

Read.Input.Data = function(setID) {
  sourceData<-list()
  sheetnames<-excel_sheets(path= here(setID, paste0("inputdata_", setID, ".xlsx")))
  
  for (i in sheetnames){
    data <- read_excel(here(setID, paste0("inputdata_", setID, ".xlsx")), sheet = i) 
    sourceData[[i]]<-data
    data <- NULL
  }
  
  return(sourceData)
}


#Input Converter for Scenario Parameters
#Processes the Scenario tab of the inputdata spreadsheet into correct formats

Convert.Scenario.Params = function(onescenario) {
  oneparamlist = as.list(onescenario)
  oneparamlist = lapply(oneparamlist, param.Convert)
  return(oneparamlist)
}

param.Convert = function(myentry) {
  if (class(myentry) == "character") {
    myparam = myentry %>%
      strsplit(split = ",") %>%
      unlist() %>%
      type.convert(as.is = TRUE)
  } else {
    myparam = myentry
  }
  
  #If there is more than one short dash in the string, format as date.
  if (class(myparam) == "character") {
    if ( length(gregexpr("-", myparam)[[1]])>1 ) {
      myparam = as.Date(myparam)
    }
  }
  
  return(myparam)
}

## Function Model.All.Scenarios() ##
#Run all scenarios in the spreadsheet and save by filegroup
#Filegroups help break up large outputs so the laptop doesn't freeze

Model.All.Scenarios = function(sourceData) {

  for (fgpname in unique(sourceData$Scenario$filegroup)) {
    
    onefgp = filter(sourceData$Scenario, filegroup == fgpname)
    
    #Convert sourceData Scenario tab to list of parameter lists
    allparams = apply(onefgp, 1, Convert.Scenario.Params)
    #Name each scenario using the id column
    names(allparams) = unlist(lapply(allparams, function(x) x$id))
    
    #Create a file name and location for the group
    myfilename = here(setID, paste0("RawData_", setID, "_", fgpname, ".RData"))
    
    #If the filegroup raw data file exists, load it
    if (file.exists(myfilename)) {
      load(myfilename)
    } else {
      #Otherwise create empty list objects to fill
      allMCs = list()
      results = list()
    }
    
    #Run each scenario if it doesn't already exist in the data file
    for (myid in names(allparams)) {
      
      if ( (myid %in% names(results)) == FALSE) {
        print(paste0("Starting computations for scenario id ", myid, 
                     " at ", format(Sys.time(), "%H:%M:%S") ))
        #Generate MC variables for each scenario
        allMCs[[myid]] = monte_carlo_vars( param_list = allparams[[myid]], 
                                           sourceData = sourceData )
        #Add background dates to each scenario parameter set
        allparams[[myid]]$bkgd_dates = allMCs[[myid]]$bkgd_dates
        
        #Run model for each scenario
        results[[myid]] = pmap(allMCs[[myid]]$MC_vars, 
                            calc_Serum_PFAS, 
                            parameter_list = allparams[[myid]], 
                            source_data = sourceData,
                            .progress = TRUE)
        
        #Save the output and the associated files after each scenario
        save(allparams, allMCs, results, 
             file = myfilename )
        
      }
    }
    
  }
  
}

## Function Extract.Serum.at.Dates() ##
##Collect serum at the dates where downstream analysis is needed into one object

Extract.Serum.at.Dates = function(setID, sourceData) {
  
  serumdatelist = list()
  pjinfo = sourceData$ProjectionInfo
  fdinfo = sourceData$FitDistInfo
  evalinfo = sourceData$EvalInfo
  if (length(pjinfo$ScenIDs) > 0) {
  idslist = strsplit(pjinfo$ScenIDs, ",") 
  } else {
    idslist = list()
  }
  fgps = unique(sourceData$Scenario$filegroup)
  
  #For each result, get name, check which dates associate with which names, 
  #then pull those dates
  
  for (onefgp in fgps) {
    print(paste0("Starting processing for ", onefgp, " at ", format(Sys.time(), "%H:%M:%S")))
    
    #Load the data for the filegroup. We use multiple groups because a 
    #laptop may struggle with a very large file
    load(here(setID, paste0("RawData_", setID, "_", onefgp, ".RData")))
    
    onefgserum = data.frame()
    for (j in 1:length(results)) {
      
      #Which dates are needed for each scenario?
      scenID = names(results)[j]

      projdates = unique(pjinfo$Date[unlist(lapply(idslist, function(x) {scenID %in% x}))])
      fitdates = unique(fdinfo$Date[which(fdinfo$ScenID == scenID)])
      evaldates = unique(evalinfo$Date[which(evalinfo$ScenID == scenID)])
      mydates = sort(unique(c(projdates, fitdates, evaldates)))
      
      for (k in 1:length(mydates)) {
        mycolname = paste0(names(results)[j], "_", mydates[k])
        
        for (i in 1:length(results[[j]])) {
          oneiter = results[[j]][[i]]
          onefgserum[i,mycolname] = oneiter$Total_Serum[which(oneiter$Date == mydates[k])]
          
        }
      }
    }
    
    serumdatelist[[scenID]] = onefgserum
    
  }
  
  iternums = unlist(lapply(serumdatelist, nrow))
  
  if ((var(iternums)==0)| (length(iternums)==1 ) ){
    serumatdates = bind_cols(serumdatelist)
    save(serumatdates, file = here(setID, paste0("SerumatDates_", setID, ".RData")) )
  } else {
    print("Iteration numbers unequal. Data is saved as list object, but will need to become dataframe for next steps.")
    save(serumdatelist, file = here(setID, paste0("SerumDateList_", setID, ".RData")) )
  }
  
}


## Function Evaluate.Against.Observed.Sera() ##
#Compare model against observations in EvalInfo tab

Evaluate.Against.Observed.Sera = function(setID, sourceData) {
  load(here(setID, paste0("SerumatDates_", setID, ".RData")) )
  evalinfo = sourceData$EvalInfo
  evalframe = data.frame()
  
  #Calculate serum metrics for each modeled scenario
  for (i in 1:nrow(evalinfo)) {
    #Attach date to ScenID so it matches column names in serumatdates
    scenwdate = paste0(evalinfo$ScenID[i], "_", evalinfo$Date[i])
    
    #Divide metrics entry into a vector
    metrics = strsplit(evalinfo$Metrics[i], ",")[[1]]
    
    #If the row corresponds to a modeled scenario and date in this Set
    #calculate the metric values
    if (scenwdate %in% colnames(serumatdates)) {
      evalvec = c()
      
      #Calculate each metric for the modeled scenario
      for (j in 1:length(metrics)) {
        if (metrics[j] == "arimean") {
          evalvec[j] = mean(serumatdates[,scenwdate])
        } else if (metrics[j] == "geomean") {
          evalvec[j] = exp(mean(log(serumatdates[,scenwdate])))
        } else {
          evalvec[j] = quantile(serumatdates[,scenwdate], as.numeric(metrics[j]))
        }
      }
      
    } else {
      #If it's not a modeled scenario, there should be metric values already in the spreadsheet
      evalvec = as.numeric(strsplit(evalinfo$Values[i], ",")[[1]])
    }
    
    #Put together the metrics and the values into a table
    #bind_rows will allow NAs where some rows don't have a particular metric available
    names(evalvec) = paste0("Metric_",metrics)
    evalframe = bind_rows(evalframe, evalvec)
    
  }
  
  #Attach the split-up information onto the original evaluation table
  evaltable = evalinfo %>%
    dplyr::select(-c(Metrics, Values)) %>%
    bind_cols(evalframe)
  
  write.csv(evaltable, file = here(setID, paste0("EvalTable_", setID, ".csv")))
  save(evaltable, file = here(setID, paste0("EvalTable_", setID, ".RData")) )
  
}


## Function Fit.Serum.Distributions() ##
#Fit logspline densities to selected scenarios at a date

Fit.Serum.Distributions = function(setID, sourceData) {
  load(here(setID, paste0("SerumatDates_", setID, ".RData")) )
  fdinfo = sourceData$FitDistInfo
  
  if ( !all((unique(fdinfo$FitType) %in% c("logspline", "none")))) {
    print("Some fit types cannot be processed. Check FitType column")
  }
  
  #Fit logspline distributions
  if ("logspline" %in% unique(fdinfo$FitType)) {
    print("Fitting logsplines to serum distributions at indicated dates")
    
    lsfinfo = filter(fdinfo, FitType == "logspline")
    logsplinefits = list()
    
    for (i in 1:nrow(lsfinfo)) {
      dataID = paste0(lsfinfo$ScenID[i], "_", lsfinfo$Date[i])
      distID = lsfinfo$DistributionID[i]
      logsplinefits[[distID]] = logspline(serumatdates[,dataID], lbound = 0)
    }
    
    save(logsplinefits, file = here(setID, paste0("LogsplineFits_", setID, ".RData")))
    
  }
  
}


## Function Project.from.Serum() ##
#Project a dependent variable (e.g. water concentration, exposure duration) based on the serum data

Project.from.Serum = function(setID, sourceData) {
  
  load(here(setID, paste0("SerumatDates_", setID, ".RData")) )
  pjinfo = sourceData$ProjectionInfo
  pjoutput = pjinfo
  scensum = sourceData$ScenSummary
  idslist = strsplit(pjinfo$ScenIDs, ",")
  
  for (i in 1:nrow(pjinfo)) {
    #Get the serum information from the serumatdates table and take relevant percentile.
    idsforproj = paste0(idslist[[i]], "_", pjinfo$Date[i])
    SerumVal = apply(serumatdates[,idsforproj], 2, quantile, probs = pjinfo$Percentile[i]/100)
    
    #Use the ScenSummary tab to find the dependent variable value for each scenario
    ivname = pjinfo$IndependentVar[i]
    idpvar = scensum[match(idslist[[i]], scensum$ScenID),ivname]
    
    mydata = cbind(idsforproj, idpvar, SerumVal)
    
    if (pjinfo$RegressionType[i] == "Linear") {
      #Make the linear fit
      myformula = as.formula(paste0("SerumVal ~ ", ivname))
      myfit = lm(formula = myformula, data = mydata)
      
      #Add projection results to the projection table
      pjoutput[i,"Projection"] = (pjinfo$SerumThreshold[i] - myfit$coefficients[1]) / myfit$coefficients[2]
      pjoutput[i,"Intercept"] = myfit$coefficients[1]
      pjoutput[i,"Slope"] = myfit$coefficients[2]
      pjoutput[i,"RSquared"] = summary(myfit)$r.squared
      
      #Plot to evaluate regression quality
      
      #Make the folder for the plots if needed
      if (dir.exists(here(setID, "ResidualPlots")) == FALSE) {
        dir.create(here(setID, "ResidualPlots"))
      }
      
      #Start the image file
      png(filename = here(setID, "ResidualPlots", paste0("LinearFitCheck_", pjoutput$ProjID[i], ".png")), 
          width = 700, height = 400,
          res = 100)
      par(mfrow = c(1,2))
      
      #Plot absolute residuals
      plot(x = mydata[,ivname], 
           y = residuals(myfit), 
           xlab = ivname, 
           ylab = "Serum Residual (ug/L)" )
      title("Absolute Residuals", cex.main = 1, line = 1, font.main = 1)
      mtext(paste0("R^2 = ", format(signif(summary(myfit)$r.squared, digits = 4), nsmall = 4) ), 
            line = 1, font = 1, side = 3, adj = 1.5, cex = 0.8)
      # Plot relative residuals
      plot(x = mydata[,ivname], 
           y = 100*residuals(myfit)/mydata$SerumVal, 
           xlab = ivname, 
           ylab = "% Serum Residual" )
      title("Relative Residuals", cex.main = 1, line = 1, font.main = 1)
      #Add main title
      mtext(paste0("Linear Fit Check for ", pjoutput$ProjID[i], " in ", setID), outer = TRUE, line = -2, font = 2)
      #Finish the image file
      dev.off()
      
    } else {
      print("Regression type not recognized. Projection cannot be completed.")
    }
    #We could add other types of fits later using this format. 
    
  }
  
  write.csv(pjoutput, 
            file = here(setID, paste0("ProjectionOutput_", setID, ".csv")),
            row.names = FALSE)
  
}


## Function Prep.Projection.Plot.Data() ##
#Generate data table for dependent variable vs serum quantiles plots

Prep.Projection.Plot.Data = function(setID, sourceData, addlquants = NA) {
  
  load(here(setID, paste0("SerumatDates_", setID, ".RData")) )
  scensum = sourceData$ScenSummary
  pjinfo = sourceData$ProjectionInfo
  
  #Assemble information for each scenario ID and date for later sorting or plotting
  scenid = substr(colnames(serumatdates), 1, nchar(colnames(serumatdates))-11)
  scendates = substr(colnames(serumatdates), nchar(colnames(serumatdates))-9, nchar(colnames(serumatdates)))
  
  quantsatdates = cbind(colnames(serumatdates), scensum[match(scenid, scensum$ScenID),])
  colnames(quantsatdates)[1] = "IDwithDate"
  quantsatdates$Date = scendates
  
  
  #Look up and glue together the PJIDs associated with each IDwithDate.
  scenidbyproj = strsplit(pjinfo$ScenIDs, ",")
  
  for (i in 1:length(scenidbyproj)) {
    scenidbyproj[[i]] = paste0(scenidbyproj[[i]], "_", pjinfo$Date[i])
  }
  
  for (i in 1:nrow(quantsatdates)) {
    
    myholder = paste0(pjinfo$ProjID[sapply(scenidbyproj, function(x) {
      quantsatdates$IDwithDate[i] %in% x
    }) ] , collapse = ",")
    
    quantsatdates$ProjIDs[i] = myholder
    
  }
  
  #Combine quantiles from defaults, projection info, and user specification
  projquants = unique(pjinfo$Percentile/100)
  myquants = c(0.05, 0.25, 0.5, 0.75, 0.95)
  myquants = unique(c(myquants, projquants, addlquants/100))
  myquants = sort(myquants[!is.na(myquants)])
  
  #Calculate all the quantiles for all scenarios and dates
  quanthold = t(apply(serumatdates, 2, quantile, probs = myquants))
  colnames(quanthold) = paste0("Serum.p", gsub("%", "", colnames(quanthold)))
  
  #Attach to the main dataframe
  quantsatdates = cbind(quantsatdates, quanthold)
  
  rownames(quantsatdates) = NULL
  
  save(quantsatdates, file = here(setID, paste0("ProjPlotQuants_", setID, ".RData")))
  write.csv(quantsatdates, file = here(setID, paste0("ProjPlotQuants_", setID, ".csv")))
  
}


## Function Prep.Lifetime.Plot.Data() ##
#Calculate all the quantiles for all scenarios and dates
#as would be needed for the lifetime serum plots

#This is built as a long table. The plots need it as a long table, 
# where the quantiles are a column. This will help with the potential problem 
# if some scenarios don't have dates that other scenarios have available

Prep.Lifetime.Plot.Data = function(setID, sourceData, addlquants = NA) {
  
  pjinfo = sourceData$ProjectionInfo
  fgps = unique(sourceData$Scenario$filegroup)
  
  #Combine quantiles from defaults, projection info, and user specification
  projquants = unique(pjinfo$Percentile/100)
  myquants = c(0.05, 0.25, 0.5, 0.75, 0.95)
  myquants = unique(c(myquants, projquants, addlquants/100))
  myquants = sort(myquants[!is.na(myquants)])
  alldatesquants = data.frame()
  
  for (onefgp in fgps) {
    print(paste0("Starting processing for ", onefgp, " at ", format(Sys.time(), "%H:%M:%S")))
    
    if (onefgp %in% sourceData$AgeInfo$FileGroup) {
      print("Cannot make lifetime plot for scenario with varying ages.")
    } else {
      
      #Load the data for the filegroup. We use multiple groups because 
      #a laptop may struggle with a very large file
      load(here(setID, paste0("RawData_", setID, "_", onefgp, ".RData")))
      
      for (j in 1:length(results)) {
        
        scenID = names(results)[j]
        onescenserum = matrix(nrow = nrow(results[[j]][[1]]), ncol = length(results[[j]]))
        
        for (i in 1:length(results[[j]])) {
          oneiter = results[[j]][[i]]
          onescenserum[,i] = oneiter$Total_Serum
        }
        
        onescenquants = t(apply(onescenserum, 1, quantile, p = myquants))
        
        colnames(onescenquants) = paste0("p", gsub("%", "", colnames(onescenquants)))
        onescenquants = cbind(results[[j]][[1]]$Date, as.data.frame(onescenquants))
        colnames(onescenquants)[1] = "Date"
        
        holder = pivot_longer(onescenquants, cols = !Date, names_to = "Percentile", values_to = "Serum")
        holder$ScenID = scenID
        alldatesquants = rbind(alldatesquants, holder)
        
      }
    }
    
  }
  
  save(alldatesquants, file = here(setID, paste0("AllDatesQuants_", setID, ".RData")) )
  write.csv(alldatesquants, file = here(setID, paste0("AllDatesQuants_", setID, ".csv")))
  
}


## Function Check.Stability() ##
#Function to check stability for a whole scenario set and make a table with results

Check.Stability = function(setID, sourceData, percentileset, error) {
  
  load(here(setID, paste0("SerumatDates_", setID, ".RData")) )
  
  stabilitytable = data.frame(ScenIDwDate = colnames(serumatdates),
                              nItersUsed = apply(serumatdates, 2, function(x) {sum(!is.na(x))})) 
  
  for (percentile in percentileset) {
    
    holder2 = apply(serumatdates, 2, quantile, probs = percentile, na.rm = TRUE)
    
    stabilitytable[, paste0("Serump",percentile*100)] = holder2
    
    for (i in 1:nrow(stabilitytable)) {
      
      print(paste0("Calculating stability at ", percentile*100, "percentile for ", stabilitytable$ScenIDwDate[i]))
      
      stabilitytable[i, paste0("StableItersp",percentile*100)] = 
        Calculate.Stable.Iterations(serumatdates[,i], 
                                    ptile = percentile, 
                                    error = error,
                                    errortype = "percentage",
                                    conf = 0.999)
      
      stabilitytable[i, paste0("PercStabilityp",percentile*100)] = 
        Calculate.Percent.Stability(serumatdates[,i], 
                                    ptile = percentile, 
                                    niter = stabilitytable$nItersUsed[i], 
                                    conf = 0.999)
    }
    
    unstable = stabilitytable[stabilitytable$nItersUsed < stabilitytable[, paste0("StableItersp",percentile*100)], ]
    
    maxiters = max(stabilitytable[,paste0("StableItersp",percentile*100)])
    maxdex = which(stabilitytable[,paste0("StableItersp",percentile*100)] == maxiters)
    
    print(paste0("The highest number of iterations required to reach ", error*100, "% stability at the ", percentile*100, "th percentile is ", maxiters))
    print("The scenario(s) with date estimated to require this number of iterations is:")
    cat(paste(stabilitytable$ScenIDwDate[maxdex], collapse = "\n"))
    cat("\n")
    print(paste0("The following scenarios with dates fail ", error*100, "% stability at the ", percentile*100, "th percentile:"))
    cat(paste(unstable$ScenIDwDate, collapse = "\n"))
    cat("\n")
  }
  
  write.csv(stabilitytable, file = here(setID, paste0("StabilityCheck_", setID, ".csv")))
  
}


## Function Calculate.RMSE() ##
#Calculate log10 root mean square error

Calculate.RMSE = function(setID, sourceData, MetricGroups = list(), 
                          ChemicalGroups = list(), PopulationGroups = list()) {
  
  load(here(setID, paste0("EvalTable_", setID, ".RData")))
  
  #RMSE is only calculated for modeled data
  logsqresids = evaltable %>%
    filter(ParmType != "Obs")
  
  #Find which columns are metrics
  metriccols = grep("Metric", colnames(logsqresids)) 
  
  for (i in 1:nrow(logsqresids)) {
    #Match the metrics in each row of the model to the metrics of appropriate observation
    matchedcols = evaltable[which(as.character(logsqresids[i,"ObservationToCompare"]) == evaltable$ScenID), metriccols]
    
    #Prepare to compare minima and maxima against 0.05 and 0.95 percentiles
    if (!is.null(matchedcols$Metric_min)) {
      if (is.na(matchedcols$Metric_0.05) & (!is.na(matchedcols$Metric_min))) {
        matchedcols$Metric_0.05 = matchedcols$Metric_min
      }
    }
    
    if (!is.null(matchedcols$Metric_max)) {
      if (is.na(matchedcols$Metric_0.95) & (!is.na(matchedcols$Metric_max))) {
        matchedcols$Metric_0.95 = matchedcols$Metric_max
      }
    }
    
    #Calculate log10 square residual for each model scenario and metric against the appropriate observation
    logsqresids[i,metriccols] = (log10(matchedcols) - log10(logsqresids[i,metriccols]))^2
  }
  
  #Make table longer and remove NA residuals where comparisons weren't available
  logsqresids = pivot_longer(logsqresids, colnames(logsqresids)[metriccols], names_to = "Metric", values_to = "LogSqResidual")
  logsqresids = logsqresids[!is.na(logsqresids$LogSqResidual),]
  
  #Prepare RMSE categories
  #Groupings are important because RMSE is a relative measure
  #so most useful to ask 'did the model fit some data better than others'
  #or 'did model parameter A or B fit the data better'
  
  #Add an 'All' group to each of the groupings
  MetricGroups = c(list("All" = unique(logsqresids$Metric)), MetricGroups)
  #Also add groupings that are just each chemical on its own
  singlechems = as.list(unique(logsqresids$Chemical))
  names(singlechems) = unique(logsqresids$Chemical)
  ChemicalGroups = c(list("All" = unique(logsqresids$Chemical)), singlechems, ChemicalGroups)
  PopulationGroups = c(list("All" = unique(logsqresids$DemoGroup)), PopulationGroups)
  
  rmsetab = data.frame()
  
  #Loop through prediction types
  #This is relevant when comparing between model parameter options to find the better fit
  for (mypt in unique(logsqresids$ParmType)) {
    oneparmtype = logsqresids %>%
      filter(ParmType == mypt)
    
    #Loop through metric, chemical, and population subgroups
    #which all contain an automatically added 'All' group
    for (i in 1:length(MetricGroups)) {
      
      for (j in 1:length(ChemicalGroups)) {
        
        for (k in 1:length(PopulationGroups)) {
          subset = oneparmtype %>%
            filter(Metric %in% MetricGroups[[i]]) %>%
            filter(Chemical %in% ChemicalGroups[[j]]) %>%
            filter(DemoGroup %in% PopulationGroups[[k]])
          
          #Calculate RMSE for each subgroup
          rmserow = data.frame(ParmType = mypt,
                               PopulationGroup = names(PopulationGroups)[k],
                               ChemicalGroup = names(ChemicalGroups)[j], 
                               MetricGroup = names(MetricGroups)[i],
                               RMSE = sqrt(mean(subset$LogSqResidual, na.rm = TRUE)))
          
          rmsetab = rbind(rmsetab, rmserow)
          
        }
      }
    }
    
  }
  
  #Save the output as a .csv file
  write.csv(rmsetab, file = here(setID, paste0("RMSEtable_", setID, ".csv")))
  
  #Return the R object so output table can be made wider.
  #For our use, we like it wider by chemical group
  return(rmsetab)
  
}
