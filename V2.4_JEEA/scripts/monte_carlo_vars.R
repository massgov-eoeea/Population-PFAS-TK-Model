### monte_carlo_vars.R ###
#Generates random parameters drawn from distributions specified in sourceData
#All functions required to make the MonteCarlo variables are now in this file

## Function make_rnd() ##
#This function makes randomly distributed values based on distribution type and parameters
#It is needed to match the distribution name in the spreadsheet to the relevant inbuilt function

#You must match the case to these distribution names.
#Lognormal requires specifying if you are using geometric mean or the meanlog

make_rnd <- function(parm1, parm2, parm3 = NA, dtype = "Normal", n = 1000){
  
  if (dtype == "Normal") {
    R = rnorm(n = n, mean = parm1, sd = parm2)
  } else if (dtype == "Lognormal_geommean") {
    R = rlnorm(n = n, meanlog = log(parm1), sdlog = log(parm2))
  } else if (dtype == "Lognormal_meanlog") {
    R = rlnorm(n = n, meanlog = parm1, sdlog = parm2)
  } else if (dtype == "Weibull") {
    R = rweibull(n = n, shape = parm1, scale = parm2)
  } else if (dtype == "Gamma") {
    R = rgamma(n = n, shape = parm1, rate = parm2)  
  } else if (dtype == "None") {
    R = rep(parm1, n)
  } else if (dtype == "Uniform") {
    R = runif(n = n, min = parm1, max = parm2)  
  } else if (dtype == "Triangular") {
    R = rtri(n = n, min = parm1, max = parm2, mode = parm3)  
  } else {
    print(paste0("Cannot match distribution type ", dtype))
    print("Check function make_rnd() in monte_carlo_vars.R file")
  }
  
  return(R)
}

## Function monte_carlo_vars() ##
#Generate iterations of parameters to use in the model

monte_carlo_vars <- function(param_list, sourceData){
  require(dplyr)
  require(purrr)
  
  lookup_tab <- sourceData[["MonteCarlo"]]
  
  #Filter table to select specific PFAS chemical parameters
  lookup_tab <- lookup_tab[lookup_tab$Chemical %in% 
                             c(NA,"NA", param_list$chemical),]
  
  #Convert table to numeric type for calculations
  lookup_tab <- lookup_tab %>% mutate(Mean = as.numeric(lookup_tab$Mean),
                                      Sd = as.numeric(lookup_tab$Sd)) %>% suppressWarnings()
  # ignore the warning about NA introduced by coercion
  
  # create MC inputs matrix
  print("Creating Monte Carlo variables dataframe")
  
  MC_vars <- data.frame("iter" = 1:param_list$niter, 
      "DWIaf_Infant" = make_rnd(dtype = lookup_tab$Distribution[which(lookup_tab$Variable == "DWIaf_Infant")],
         n = param_list$niter, 
         parm1 = lookup_tab$Mean[which(lookup_tab$Variable == "DWIaf_Infant")],
         parm2 = lookup_tab$Sd[which(lookup_tab$Variable == "DWIaf_Infant")]),
      "DWIaf_ChildAdult" = make_rnd(dtype = lookup_tab$Distribution[which(lookup_tab$Variable == "DWIaf_ChildAdult")],
         n = param_list$niter, 
         parm1 = lookup_tab$Mean[which(lookup_tab$Variable == "DWIaf_ChildAdult")],
         parm2 = lookup_tab$Sd[which(lookup_tab$Variable == "DWIaf_ChildAdult")]),
      "DWIaf_Teen" = make_rnd(dtype = lookup_tab$Distribution[which(lookup_tab$Variable == "DWIaf_Teen")],
         n = param_list$niter, 
         parm1 = lookup_tab$Mean[which(lookup_tab$Variable == "DWIaf_Teen")],
         parm2 = lookup_tab$Sd[which(lookup_tab$Variable == "DWIaf_Teen")]),
      "HalfLife" = make_rnd(dtype = lookup_tab$Distribution[which(lookup_tab$Variable == "HalfLife")],
         n = param_list$niter, 
         parm1 = lookup_tab$Mean[which(lookup_tab$Variable == "HalfLife")],
         parm2 = lookup_tab$Sd[which(lookup_tab$Variable == "HalfLife")]),
      "VolDist" = make_rnd(dtype = lookup_tab$Distribution[which(lookup_tab$Variable == "VolDist")],
          n = param_list$niter, 
          parm1 = lookup_tab$Mean[which(lookup_tab$Variable == "VolDist")],
          parm2 = lookup_tab$Sd[which(lookup_tab$Variable == "VolDist")]),
      "LactXfer" = make_rnd(dtype = lookup_tab$Distribution[which(lookup_tab$Variable == "LactXfer")],
          n = param_list$niter, 
          parm1 = lookup_tab$Mean[which(lookup_tab$Variable == "LactXfer")],
          parm2 = lookup_tab$Sd[which(lookup_tab$Variable == "LactXfer")]), 
      "PlacXfer" = make_rnd(dtype = lookup_tab$Distribution[which(lookup_tab$Variable == "PlacXfer")],
         n = param_list$niter, 
         parm1 = lookup_tab$Mean[which(lookup_tab$Variable == "PlacXfer")],
         parm2 = lookup_tab$Sd[which(lookup_tab$Variable == "PlacXfer")]),
      "BirthLoss" = make_rnd(dtype = lookup_tab$Distribution[which(lookup_tab$Variable == "BirthLoss")],
          n = param_list$niter, 
          parm1 = lookup_tab$Mean[which(lookup_tab$Variable == "BirthLoss")] ,
          parm2 = lookup_tab$Sd[which(lookup_tab$Variable == "BirthLoss")]  ) ,
      "kMilk" = make_rnd(dtype = lookup_tab$Distribution[which(lookup_tab$Variable == "kMilk")],
         n = param_list$niter, 
         parm1 = lookup_tab$Mean[which(lookup_tab$Variable == "kMilk")],
         parm2 = lookup_tab$Sd[which(lookup_tab$Variable == "kMilk")]) 
      )

  #Prepare a centralLX to compare to each iteration's LactXfer
  #Set a default centralLX based on the main table
  centralLX = lookup_tab$Mean[which(lookup_tab$Variable == "LactXfer")]
  
  #But if an alternative LX gets used, we need to use that one instead. 

  #Replace relevant columns with alternative parameters if the ScenarioID appears in the AltParams tab
  altlookuptab = sourceData[["AltParams"]]
  
  if (param_list$id %in% altlookuptab$ScenarioID) {
    altvars = altlookuptab %>%
      filter(ScenarioID == param_list$id)
    
    for (i in 1:nrow(altvars)) {
      varname = altvars$Variable[i]
      MC_vars[,varname] = make_rnd(parm1 = altvars$Parm1_Mean[i],
                                      parm2 = altvars$Parm2_Sd[i], 
                                      parm3 = altvars$Parm3[i], 
                                      dtype = altvars$Distribution[i],
                                      n = nrow(MC_vars))
      
      #Store the alternative LactXfer mean if an alternative LactXfer was used
      if(varname == "LactXfer") {
        centralLX = altvars$Parm1_Mean
      }
    
    }
  
  }
  
  #Adjust kMilk by the LactXfer of each iteration relative to the central LX
  #Even if we calculate kMilk from a known milk decline over time, the variation
  #in lactational transfer should connect to the breastmilk elimination rate constant. 

  MC_vars = MC_vars %>%
    mutate(kMilk = kMilk*LactXfer/centralLX)
  #mutate() adjusts kMilk in each row of the MC_vars table 
  #using the original values in that row for columns kMilk and LactXfer
  #as well as centralLX, a pre-defined constant outside the MC_vars table. 

  #Generate starting maternal serum concentrations
  
  if (param_list$mom_serum_method == "lookup") {
    #Use maternal serum distribution from MatSerum tab
    momserum_tab <- sourceData[["MatSerum"]]
    
    momserum_tab <- momserum_tab %>%
      filter(Chemical %in% c(NA,"NA", param_list$chemical) ) %>%
      filter(ParmType == "TotalSerum") %>%
      mutate(meanlog = as.numeric(meanlog),
             sdlog = as.numeric(sdlog)) %>% 
      suppressWarnings()
    
    #Use the row in MatSerum where WaterID matches mom_water_conc 
    MC_vars[,"MatStart"] = make_rnd(dtype = momserum_tab$Distribution[which(momserum_tab$WaterID == param_list$mom_water_conc)], 
        n = param_list$niter, 
        parm1 = momserum_tab$meanlog[which(momserum_tab$WaterID == param_list$mom_water_conc)],
        parm2 = momserum_tab$sdlog[which(momserum_tab$WaterID == param_list$mom_water_conc)]  )
    #MatStart is in ug/L
    
  } else if (param_list$mom_serum_method == "steadystate") {
    #Use parameters to create a steady-state maternal distribution
    
    #Draw the random maternal background from the MatSerum tab BkgdSerum rows
    mombkgd_dist <- sourceData[["MatSerum"]]
    
    mombkgd_dist <- mombkgd_dist %>%
      filter(Chemical %in% c(NA,"NA", param_list$chemical) ) %>%
      filter(ParmType == "BkgdSerum") %>%
      mutate(meanlog = as.numeric(meanlog),
             sdlog = as.numeric(sdlog)) %>% 
      suppressWarnings()
    
    MC_vars[,"MatBkgd"] = make_rnd(dtype = mombkgd_dist$Distribution, 
                                    n = param_list$niter, 
                                    parm1 = mombkgd_dist$meanlog,
                                    parm2 = mombkgd_dist$sdlog )
    
    #Incorporate the relative source contribution factor of 0.8
    MC_vars[,"MatBkgd"] <- MC_vars[,"MatBkgd"]*0.8   
    
    #Using EFH Table 3-63 mean water intake for women of childbearing age = 15.6 mL/(kg*day)
    MC_vars[,"MatStart"] = (( (param_list$mom_water_conc/1000) * 0.0156 * MC_vars$DWIaf_ChildAdult * MC_vars$HalfLife *365.25 ) /( MC_vars$VolDist * log(2))) + MC_vars$MatBkgd
    
    #To save space, we don't save MatBkgd in the output
    MC_vars[,"MatBkgd"] = NULL
    
  } else if (param_list$mom_serum_method == "logspline") {
    #Look up maternal serum distributions from a set of logspline fits
    
    #They should be appended the sourceData list in a named list themselves
    momlogsplines = sourceData[["MatLogsplines"]]
    
    #where name corresponds to the mom_water_conc used
    MC_vars[,"MatStart"] = rlogspline(n = param_list$niter, 
               fit = momlogsplines[[as.character(param_list$mom_water_conc)]])
  }
  
  #Generate water concentrations
  
  if (is.null(param_list$water_concs_sds)) {
    #If there isn't a column for water conc standard deviation, prepare to use constant water conc
    
    water_concs_to_nest <- data.frame(parm1 = as.numeric(param_list$water_concs), 
                    parm2 = NA,
                    dtype = "None", 
                    n = rep(param_list$niter, length(param_list$water_concs)), 
                    stringsAsFactors = FALSE)
    
  } else {
    #Prepare water concentration distributions
    water_concs_to_nest <- data.frame(parm1 = as.numeric(param_list$water_concs), 
                    parm2 = as.numeric(param_list$water_concs_sds),
                    dtype = as.character(param_list$water_concs_dists), 
                    n = rep(param_list$niter, length(param_list$water_concs)), 
                    stringsAsFactors = FALSE)
  }
  
  #Draw water concentrations from pre-specified distribution
  water_concs_to_nest<- pmap(water_concs_to_nest,make_rnd) %>% as.data.frame()
  
  #Water conc cannot be negative. Important for GustavusUS evaluation
  water_concs_to_nest[water_concs_to_nest<0] = 0  
  
  #Finish setting up water conc holder
  colnames(water_concs_to_nest) <- paste0("t", 1:ncol(water_concs_to_nest))
  water_concs_to_nest <- water_concs_to_nest %>% 
    mutate("iter" = 1:param_list$niter)  
  
  
  #Create age distribution if specified
  
  agetab = sourceData[["AgeInfo"]]
  
  if (param_list$filegroup %in% agetab$FileGroup) {
  #Prepare distribution table with the number of iterations for each age group
    scenages = agetab %>%
      filter(FileGroup == param_list$filegroup) %>%
      mutate(Fraction = as.numeric(Fraction), 
             Param1 = as.numeric(Param1),
             Param2 = as.numeric(Param2),
             Param3 = as.numeric(Param3)) %>%
      mutate(Niter = ifelse( is.na(Fraction), 
                             NA,
                             round(Fraction*param_list$niter, digits = 0)) )
    
    #We use the NA in the Fraction column to denote the category that will
    #take the remainder of the other categories, to avoid problems with rounding
    remainderamt = param_list$niter - sum(scenages$Niter, na.rm = TRUE)
    if (remainderamt > 0) {
      remainderloc = which(is.na(scenages$Fraction))
      scenages$Niter[remainderloc] = remainderamt
    }
    #This ensures that the total number of start dates generated is exactly equal to 
    #the total number of iterations we will run
    
    #Draw ages for each age group in the scenario
    MCages = numeric()
    
    for (i in 1:nrow(scenages)) {
      holder = floor(make_rnd(parm1 = scenages$Param1[i], 
                              parm2 = scenages$Param2[i], 
                              parm3 = scenages$Param3[i],
                              dtype = scenages$Distribution[i], 
                              n = scenages$Niter[i]
      ))
      MCages = c(MCages, holder)
      
    }
    
    #Convert the drawn ages into start dates (aka birthdays)
    styears = year(param_list$ed)-MCages
    #All iterations share a birthday as dictated by the original start date,
    #but are given different birth years depending on age
    MCstartdates = param_list$sd
    #Randomize order of start_dates so the first entries are not youngest
    year(MCstartdates) = sample(styears)
    
    #Store the varied start dates as a column in the MC_vars dataframe
    MC_vars$StartDate = MCstartdates
    
  } 

  
  #Generate background serum concentrations
  
  bkgd_tab <- sourceData[["BkgdSerum"]]
  
  #Filter source by chemical and convert date column
  bkgd_tab = bkgd_tab %>%
    filter(Chemical == param_list$chemical) %>%
    mutate(StartDate = as.Date(StartDate)) %>%
    as.data.frame()
  
  #Select receptor gender
  if (param_list$gender == "Female") {
    groupID = c("BkgdFC", "BkgdFA")
  } else if (param_list$gender == "Male") {
    groupID = c("BkgdMC", "BkgdMA")
  } else if (param_list$gender == "Unknown") {
    groupID = c("BkgdUC", "BkgdUA") 
  }
  
  #Find the relevant date intervals and store in bkgd_dates
  
  #Check if receptor is under age 12
  bday12yo = add_with_rollback(param_list$sd, years(12))  #Adult is defined as being over age 12
  
  if (param_list$ed < bday12yo) {
    childrows = findInterval(param_list$sd, bkgd_tab$StartDate) : findInterval(param_list$ed, bkgd_tab$StartDate)
    adultrows = c()
    bkgd_dates = bkgd_tab[childrows, "StartDate"]
  } else {
    childrows = findInterval(param_list$sd, bkgd_tab$StartDate) : findInterval(bday12yo, bkgd_tab$StartDate)
    adultrows = findInterval(bday12yo, bkgd_tab$StartDate) : findInterval(param_list$ed, bkgd_tab$StartDate)
    #The date assignment should work even if the receptor turns 12 on a date where the background changes
    bkgd_dates = c(bkgd_tab[childrows, "StartDate"], bday12yo, bkgd_tab[adultrows[-1], "StartDate"])
  }

  #Find the relevant background distributions and generate random serum background values
  bkgd_to_nest <- data.frame(
     parm1 = c(bkgd_tab[childrows, paste0(groupID[1], "_Mean")], 
               bkgd_tab[adultrows, paste0(groupID[2], "_Mean")])  , 
     parm2 = c(bkgd_tab[childrows, paste0(groupID[1], "_Sd")], 
               bkgd_tab[adultrows, paste0(groupID[2], "_Sd")])  ,
     dtype = c(bkgd_tab[childrows, "Distribution"], 
               bkgd_tab[adultrows, "Distribution"]) , 
     n = param_list$niter , 
            stringsAsFactors = FALSE)
  
  bkgd_to_nest<- pmap(bkgd_to_nest,make_rnd) %>% as.data.frame()
  colnames(bkgd_to_nest) <- paste0("t", 1:ncol(bkgd_to_nest))
  #Incorporate the relative source contribution
  #Assumed 0.8 of unexposed population's PFAS come from non-drinking-water sources
  bkgd_to_nest <- bkgd_to_nest*0.8   
  bkgd_to_nest$iter <- 1:param_list$niter

  #Join the nested concentrations into the main dataframe

  MC_vars <- MC_vars %>% 
    nest_join(as.data.frame(water_concs_to_nest), by = "iter", name = "WaterConcs") %>%
    nest_join(bkgd_to_nest, by = "iter", name = "BkgdConcs")
  
  #If using varied ages, background can only change with age, not also with time
  #We replace bkgd_dates with the keyphrase "UseAgeInstead"
  if (is.null(MC_vars$StartDate) == FALSE) {
    bkgd_dates = "UseAgeInstead"
    
    #The "UseAgeInstead" bkgd_dates option uses at most two background values - one child, one adult
    #So if there are more background values for an age-varying scenario, print a warning.
    if (ncol(bkgd_to_nest) > 3) {
      print("More than two background values are specified for a scenario with age variation")
      print("Serum calculations may produce unexpected outputs. Reset scenario to use constant background over time")
    }
  }
  
  #Return the bkgd_dates as well as MC_vars
  MClist = list(MC_vars, bkgd_dates)
  names(MClist) = c("MC_vars", "bkgd_dates")
  return(MClist)
}

