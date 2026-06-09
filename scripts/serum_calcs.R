### serum_calcs.R ###
#Calculates serum PFAS over time for one Monte Carlo iteration (parameter set)

## Function calc_Serum_PFAS() ##
#The main function putting all subsequent functions together
#Written to be compatible with pmap() across MC_vars dataframe

#For pmap to work, all column names of the MC_vars dataframe must appear as variables. 
#Which is why iter appears as an input here despite not being used. 

calc_Serum_PFAS <- function(
                parameter_list, source_data, DWIaf_Infant, DWIaf_ChildAdult, 
                DWIaf_Teen, HalfLife, VolDist, LactXfer, PlacXfer, BirthLoss,
                kMilk, MatStart, WaterConcs, BkgdConcs, StartDate = NA, iter) {
  
  if (is.na(StartDate)) {
    #When a constant age is used, the first date is taken directly from parameter_list
    result_frame <- create_frame(sd = parameter_list$sd, 
                         ed = parameter_list$ed, 
                         timezoom_start = parameter_list$timezoom_start, 
                         timezoom_end = parameter_list$timezoom_end, 
                         timezoom_interval = parameter_list$timezoom_interval, 
                         water_concs_dates = parameter_list$water_concs_dates, 
                         saved_dates = parameter_list$saved_dates)
    
  } else {
    #If varied ages are used, StartDate will be a column in the MC_vars
    #The dates evaluated in each iteration will be created using that iteration's StartDate
    
    result_frame <- create_frame(sd = StartDate, 
                         ed = parameter_list$ed, 
                         timezoom_start = parameter_list$timezoom_start, 
                         timezoom_end = parameter_list$timezoom_end, 
                         timezoom_interval = parameter_list$timezoom_interval, 
                         water_concs_dates = parameter_list$water_concs_dates, 
                         saved_dates = parameter_list$saved_dates)
  }
  
  result_frame <- pop_Liquid_Intake(data = result_frame, 
                                    source = source_data, 
                                    AF_Infant = DWIaf_Infant, 
                                    AF_ChildAdult = DWIaf_ChildAdult, 
                                    AF_Teen = DWIaf_Teen)
  
  #Watch out for list handling with WaterConcs - see notes
  result_frame <- pop_Liquid_Conc(data = result_frame, 
                          water_concs = WaterConcs, 
                          water_concs_dates = parameter_list$water_concs_dates, 
                          bfmonths = parameter_list$fed_mo, 
                          momserum = MatStart, 
                          lactxfer = LactXfer, 
                          kmilk = kMilk, 
                          birthloss = BirthLoss, 
                          placxfer = PlacXfer)
  
  #Watch out for list handling with BkgdConcs - see notes
  result_frame <- pop_Bkgd_Serum(data = result_frame, 
                                 bkgd_concs = BkgdConcs, 
                                 bkgd_dates = parameter_list$bkgd_dates )

  result_frame <- calc_Total_Serum(data = result_frame,
                                   HL = HalfLife,
                                   VD = VolDist,
                                   MatStart = MatStart, 
                                   PlacXfer = PlacXfer )
  
  if (is.character(parameter_list$saved_dates)) {
    return(result_frame)
  } else if ( is.null(parameter_list$saved_dates) ) {
    return(result_frame)
  } else if ( any(parameter_list$saved_dates %in% result_frame$Date) ) {
    
    #Any dates that match return serum. Duplicate dates are removed this way too
    result_vec = filter(result_frame, Date %in% parameter_list$saved_dates)$Total_Serum
    names(result_vec) = filter(result_frame, Date %in% parameter_list$saved_dates)$Date
    
    return(result_vec)

  } else {
    return(result_frame)
  }
  
  
}

## Function create_frame() ##
#Set up dates for evaluation

create_frame <- function(sd, ed, timezoom_start, timezoom_end, 
                         timezoom_interval, water_concs_dates, saved_dates){

  require(dplyr)
  require(lubridate)
  
  #Check that the receptor is old enough to need adultdates
  #Generate monthly dates for infancy and yearly dates beyond that
  if (ed-sd > 365){
    infantdates = seq(sd, add_with_rollback(sd, years(1)), by = "1 month")
    adultdates = c(seq(add_with_rollback(sd, years(1)), ed, by = "1 year"), ed)
  } else {
    infantdates = c(seq(sd, ed, by = "1 month"), ed)
    adultdates = NA
  }
  #add_with_rollback() is a special function for handling dates so months come out OK. 
  
  
  #Generate dates for the zoomed interval, if specified
  if(is.na(timezoom_start) | is.na(timezoom_end)) {
    zoomdates = NA
  } else {
    zoomdates = c(seq(timezoom_start, timezoom_end, by = timezoom_interval), timezoom_end)
  }
  
  #Add dates where the water changes concentration, as long as they're within the modeled time
  waterdates = water_concs_dates[between(water_concs_dates, sd, ed)]
  
  #NAs are removed automatically when the dates are joined together
  if (class(saved_dates) == "Date") {
    #Add the saved_dates to all the other dates, if specified
    alldates = unique(sort(c(infantdates, adultdates, zoomdates, waterdates, saved_dates)))
  }
  else { alldates = unique(sort(c(infantdates, adultdates, zoomdates, waterdates))) }
  
  data <- data.frame("id" = 1:length(alldates),
                     "Date" = alldates)
  
  #Generate a rounded-down age in months for first year of life, 
  #to enable liquid intake matching
  
  agemonths = interval(data$Date[1], data$Date) %>%
    time_length(unit = "months") %>%
    floor
  
  agemonths = round(agemonths/12, digits = 3)
  
  #Generate a rounded-down age in years beyond first year to match liquid intake
  
  ageyears = floor(time_length(interval(data$Date[1], data$Date), unit = "years"))
  
  data = data %>%
    mutate("Age_Years" = ageyears) %>%
    mutate("Age_Years" = ifelse(Age_Years < 1, agemonths, ageyears))
  
  return(data)
  
}


## FUnction pop_Liquid_Intake() ##
#Fill in daily liquid intake 

pop_Liquid_Intake <- function(data, source, AF_Infant = 1, AF_ChildAdult = 1, AF_Teen = 1){
  
  # Load liquid intake
  water_intake <- source[["DWI"]] %>% 
    as.data.frame() %>% 
    filter( Age_Years <= 100) # 101 was added for extrapolating
  
  #Match the DWI using the Age_Years column
  data = left_join(data, water_intake[,c("Age_Years", "DWI_mL_kg_d")], by = "Age_Years")
  
  #Multiply each liquid intake by the correct adjustment factor based on receptor age.
  data = data %>%
    mutate("Liquid_Intake" = case_when( 
      between(Age_Years, 0, 0.99) ~ DWI_mL_kg_d*AF_Infant,
      between(Age_Years, 1, 10.99) ~ DWI_mL_kg_d*AF_ChildAdult,
      between(Age_Years, 11, 20.99) ~ DWI_mL_kg_d*AF_Teen,
      between(Age_Years, 21, 101) ~ DWI_mL_kg_d*AF_ChildAdult ), 
    DWI_mL_kg_d = NULL)
  
  #The gaps in this between() call are OK because age will always be rounded
  #down to the nearest month. 

  return(data)
}


## Function pop_Liquid_Conc() ##
#Fill in PFAS concentration found in the liquid the receptor is consuming 
#at a particular date

pop_Liquid_Conc <- function(data, water_concs, water_concs_dates, bfmonths, 
                            momserum, lactxfer, kmilk, birthloss, placxfer) {
  
  water_concs_dates = as.Date(water_concs_dates)
  water_concs = unlist(water_concs)   #This is needed because of how pmap works
  
  #Assign liquid concentration from the drinking water concentration to start
  data = data %>%
    mutate("Liquid_Conc" = water_concs[findInterval(Date,water_concs_dates)])
  
  #Calculate breast milk concentration 
  
  if (bfmonths > 0) {
    #What dates are in within the breast feeding timeframe?
    milkconcs = data.frame("Date" = data$Date[data$Age_Years< round(bfmonths/12, digits = 3)])
    
    #Starting concentration of breast milk
    #Third term accounts for fraction of PFAS lost during birth, if used 
    milkconcs$BMconc[1] = momserum*1000 * lactxfer * (1-(birthloss*placxfer))
    
    #Liquid conc must be in ng/L but momserum in ug/L so multiply by 1000
    
    #Conc in milk decays exponentially based on how many days have passed
    for (i in 2:nrow(milkconcs)) {
      daynum = as.numeric(milkconcs$Date[i] - milkconcs$Date[i-1])
      milkconcs$BMconc[i] = milkconcs$BMconc[i-1] * exp(-kmilk*daynum)
    }
    
    #Use milk concentrations as the liquid concentration when breastfeeding
    data$Liquid_Conc[1:nrow(milkconcs)] = milkconcs$BMconc
  }
  
  return(data)
}

## Function pop_Bkgd_Serum() ##
#Fill in the background serum concentration for the receptor by date

pop_Bkgd_Serum <- function(data, bkgd_concs, bkgd_dates) {
  
  #If we are using variable ages... 
  if (as.character(bkgd_dates[1]) == "UseAgeInstead") {
    bkgd_concs = unlist(bkgd_concs)   #This is needed because of how pmap works
    
    #Assign background concentration based on being under or over age 12
    data = data %>%
      mutate("Bkgd_Serum" = ifelse(Age_Years < 12, bkgd_concs[1], bkgd_concs[2])) %>%
      #Background ramps up linearly from 0 for first 3 years of life
      mutate("Bkgd_Serum" = ifelse(Age_Years < 3, Bkgd_Serum*((Date-data$Date[1])/(3*365.25)), Bkgd_Serum) )
    
  } else {
    bkgd_dates = as.Date(bkgd_dates)
    bkgd_concs = unlist(bkgd_concs)   #This is needed because of how pmap works
    
    #Assign background concentration based on the date intervals
    data = data %>%
      mutate("Bkgd_Serum" = bkgd_concs[findInterval(Date,bkgd_dates)]) %>%
      #Background ramps up linearly from 0 for first 3 years of life
      mutate("Bkgd_Serum" = ifelse(Age_Years < 3, Bkgd_Serum*((Date-data$Date[1])/(3*365.25)), Bkgd_Serum) )
  }
  
  return(data)
}
#The background can only change as smoothly as the time intervals allow.


## Function nextserum() ##
#Helper function for calculating serum after each time interval

nextserum <- function(k, VD, stdate, endate, liqin, liqconc, bkgdold, bkgdnew, stser){
  
  t = as.numeric(endate-stdate)
  
  #The portion of the starting serum that is from the old background should not
  #be subject to elimination. This prevents negative term2 if bkgd rises
  term2 = stser-bkgdold

  #This is Css in Chiu 2022
  term3 = (liqin*liqconc)/(k*VD*1000*1000)
  #liqconc is in ng/L, but serums are in ug/L. Divide by 1000.
  #liqin is in mL/(kg*day), but liqconc is in ng/L, not ng/mL. Divide by 1000 again. 
  
  endser = bkgdnew + term2*exp(-k*t) + term3*(1-exp(-k*t))
  #This is the integrated rate law for serum after t days pass, if the 
  #parameters such as intake and water conc remain constant during that time.
  
  #Serum output units will be in ug/L
  
  return(endser)
  
}


## Function calc_Total_Serum() ##
#Calculate total serum using the data table

calc_Total_Serum <- function(data, HL, VD, MatStart, PlacXfer) {
  
  #Calculate rate constant k in days
  
  k = log(2)/(HL*365.25)
  
  data$Total_Serum = NA
  
  #Calculate serum PFAS that receptor is born with
  data$Total_Serum[1] = MatStart*PlacXfer
  
  #Calculate serum PFAS using the parameters from each interval. 
  
  data$Total_Serum[2] = nextserum(k = k,
                                    VD = VD, 
                                    stdate = data$Date[1],
                                    endate = data$Date[2],
                                    liqin = data$Liquid_Intake[1],
                                    liqconc = data$Liquid_Conc[1], 
                                    bkgdold = 0,
                                    bkgdnew = data$Bkgd_Serum[1],
                                    stser = data$Total_Serum[1])
  
  for (i in c(2: (nrow(data)-1) )) {
    data$Total_Serum[i+1] = nextserum(k = k,
                                      VD = VD, 
                                      stdate = data$Date[i],
                                      endate = data$Date[i+1],
                                      liqin = data$Liquid_Intake[i],
                                      liqconc = data$Liquid_Conc[i], 
                                      bkgdold = data$Bkgd_Serum[i-1],
                                      bkgdnew = data$Bkgd_Serum[i],
                                      stser = data$Total_Serum[i])
  }
  
  return(data)
  
}
