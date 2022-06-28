library(MASS)
library(tidyverse)
library(broom)
library(boot)
library(parallel)


#### data simulation functions ####
simulate_data <- function(n_cases = 1000, inflection_point = 8, max_data_duration = 24, controls = F, control_ratio = 5, n_practices = 100, age_range = 18:89, patient_var = 0.02, prac_var = 0.005, inflection_coefficient = 0) {
  # n_cases is number of cases to generate data for
  # inflection_point is number of time units pre-diagnosis that inflection point should occur
  # max_data_duration is max number of time units that a patient should have data for
  # controls is whether to produce "matched" controls
  # control_ratio is the number of controls to produce for each case if controls is True
  
  # first make sure there are no problems with entered values
  bad_time_point <- function(time_point) {
    #time points should be positive whole numbers
    bad <- F
    if (time_point < 0) {
      bad <- T
    } else if (time_point%%1 != 0) {
      bad <- T
    }
    return(bad)
  }
  
  if (bad_time_point(max_data_duration) | max_data_duration == 0) {
    stop(paste("The selected maximum data duration (", max_data_duration, "time units) is not possible."))
  } else if (bad_time_point(inflection_point)) {
    stop(paste("The selected inflection point (", inflection_point, "time units) is not possible."))
  } else if (inflection_point == 0) {
    # this may be intentional so will allow function to carry on 
    message(paste("************ WARNING ************"))
    message(paste("The inflection point (", inflection_point, "time units) occurs at diagnosis."))
    message(paste("True inflection point will not be detectable in the generated data."))
    message(paste("*********************************"))
  } else if (inflection_point >= max_data_duration) {
    # this may be intentional so will allow function to carry on
    message(paste("************ WARNING ************"))
    message(paste("The inflection point (", inflection_point, "time units) occurs before the maximum data duration allowed (", max_data_duration, "time units)."))
    message(paste("True inflection point will not be detectable in the generated data."))
    message(paste("*********************************"))
  }
  
  if ((n_cases <= 0)|(n_cases%%1 != 0)) {
    stop(paste("The selected number of cases (", n_cases, ") is not possible."))
  } else if (!rapportools::is.boolean(controls)) {
    stop(paste("The selected value of controls (", controls, ") is not possible."))
  } else if (controls) {
    if ((control_ratio%%1 != 0) | (control_ratio <= 0)) {
      stop(paste("Controls have been requested but the selected control_ratio (", control_ratio, ") is not possible."))
    }
  }
  
  # now can start building data
  
  #regardless of controls, need to generate case data
  #say latest diag_year is 2001 because require at least one year of follow up
  patients <- tibble(patid = seq(1, n_cases), 
                     pracid = rdunif(n_cases, n_practices), 
                     case_control = "case", 
                     sex = wakefield::sex(n = n_cases), 
                     age = wakefield ::age(n = n_cases, x = age_range),
                     diag_year = 2019 - rdunif(n_cases, 0, 18), 
                     diag_month = sample(1:12, n_cases, replace = T), 
                     max_follow_up = (diag_year - 2000)*12 + diag_month
  ) %>%
    mutate(female = case_when(sex == "Female" ~ 1, T ~ 0))
  
  if (controls) {
    patients$matched_case <- patients$patid
    
    #now create controls
    patients <- patients %>% 
      slice(rep(1:n(), control_ratio + 1)) %>%
      mutate(patid = seq(1, n_cases*(control_ratio + 1)),
             case_control = c(rep("case", n_cases), rep("control", n_cases*control_ratio)))
  } else {
    #set control_ratio = 0 for simplicity
    control_ratio <- 0
  }
  
  #now randomly assign each patient a follow-up length between 12 (min-follow-up length unless max < 12) and their max
  #but if longer than max data duration then cut off
  rfu <- function(m) sample(12:m, 1)
  pairmin <- function(k) min(k, max_data_duration)
  patients <- patients %>%
    mutate(follow_up = purrr::pmap_dbl(list(max_follow_up), rfu)) %>%
    mutate(follow_up = purrr::pmap_dbl(list(follow_up), pairmin)) %>%
    dplyr::select(-max_follow_up)
  
  #add patient-level and practice-level effects
  patients$pat_cluster <- rnorm((control_ratio + 1)*n_cases, 0, patient_var)
  
  prac_effect <- function(i) {
    set.seed(i)
    rnorm(1, 0, prac_var)
  }
  patients <- patients %>%
    mutate(prac_cluster = purrr::pmap_dbl(list(pracid), prac_effect),
           birth_month = sample(1:12, (control_ratio + 1)*n_cases, replace = T))
  
  #now we have a table of patients
  # need to generate a row for each month of follow-up and produce consultation counts
  if (controls) {
    consultations <- tibble(patid = rep(patients$patid, patients$follow_up), 
                            pracid = rep(patients$pracid, patients$follow_up),
                            pat_cluster = rep(patients$pat_cluster, patients$follow_up), 
                            prac_cluster = rep(patients$prac_cluster, patients$follow_up), 
                            female = rep(patients$female, patients$follow_up), 
                            age = rep(patients$age, patients$follow_up), 
                            birth_month = rep(patients$birth_month, patients$follow_up),
                            case_control = rep(patients$case_control, patients$follow_up), 
                            matched_case = rep(patients$matched_case, patients$follow_up), 
                            months_pre_diag = sequence(patients$follow_up), 
                            diag_year = rep(patients$diag_year, patients$follow_up), 
                            diag_month = rep(patients$diag_month, patients$follow_up))
  } else {
    consultations <- tibble(patid = rep(patients$patid, patients$follow_up), 
                            pracid = rep(patients$pracid, patients$follow_up),
                            pat_cluster = rep(patients$pat_cluster, patients$follow_up), 
                            prac_cluster = rep(patients$prac_cluster, patients$follow_up), 
                            female = rep(patients$female, patients$follow_up), 
                            age = rep(patients$age, patients$follow_up),  
                            birth_month = rep(patients$birth_month, patients$follow_up),
                            months_pre_diag = sequence(patients$follow_up), 
                            diag_year = rep(patients$diag_year, patients$follow_up), 
                            diag_month = rep(patients$diag_month, patients$follow_up))
  }
  
  get_current_date <- function(diag_month, diag_year, n_months, return_val = NA) {
    current_month <- diag_month
    current_year <- diag_year
    for (i in 1:n_months) {
      current_month <- current_month - 1
      if (current_month == 0) {
        current_month <- 12
        current_year <- current_year - 1
      }
    }
    
    if (return_val == "year") {
      return(current_year)
    } else if (return_val == "month") {
      return(current_month)
    } else {
      return(list("year" = current_year, "month" = current_month))
    }
  }
  
  consultations <- consultations %>%
    mutate(current_year = purrr::pmap_dbl(list(diag_month, diag_year, months_pre_diag, "year"), get_current_date), 
           current_month = purrr::pmap_dbl(list(diag_month, diag_year, months_pre_diag, "month"), get_current_date))
  
  if (controls) {
    consultations <- consultations %>% 
      mutate(alt_month = case_when((months_pre_diag <= inflection_point) & (case_control == "case") ~ inflection_point - months_pre_diag, TRUE ~ 0))
  } else {
    consultations <- consultations %>% 
      mutate(alt_month = case_when((months_pre_diag <= inflection_point) ~ inflection_point - months_pre_diag, TRUE ~ 0))
  }
  
  #need to get current age
  get_current_age <- function(current_month, current_year, diag_year, birth_month, age) {
    if (current_month < birth_month) {
      return(age - (diag_year - current_year) - 1) 
    }
    else {
      return(age - (diag_year - current_year))
    }
  }
  
  consultations <- consultations %>%
    mutate(age = purrr::pmap_dbl(list(current_month, current_year, diag_year, birth_month, age), get_current_age))
  
  #generate a person-month variable for consistency
  consultations$py <- 1
  
  #don't want size of change to be proportional to length of DW 
  #if inflection_coefficient isn't specified, and there should be an inflect
  #i.e. for earlier inflection points need to shrink the alt_month coefficient
  if ((inflection_coefficient == 0) & (inflection_point)) {
    inflection_coefficient <- 2/inflection_point
  } 
  
  # simulate some data based on round-ish numbers
  # follow a poisson distribution based around the patient-level mean defined by exp of variables
  consultations$n_consultations <- rpois(nrow(consultations), exp(
    -1
    -1e-3*0.05*consultations$current_year   #secular year-on-year increase
    # - 0.005*consultations$months_pre_diag   
    + inflection_coefficient*consultations$alt_month
    + rnorm(nrow(consultations), 0, 0.05)
    + consultations$pat_cluster
    + consultations$prac_cluster
    + 0.2*consultations$female
    + 0.005*consultations$age  # increase based on age
    + log(consultations$py)
  ))
  
  if (controls) {
    return_data <- consultations %>% dplyr::select(patid, pracid, female, age, case_control, matched_case, current_month, current_year, months_pre_diag, n_consultations)
  } else {
    return_data <- consultations %>% dplyr::select(patid, pracid, female, age, current_month, current_year, months_pre_diag, n_consultations)
  }
  
  return(return_data)
}


#### GA inf point method ####

#2 alternative implementations
GA_basic_method <- function(data, controls = F) {
  # as in https://bjgp.org/content/72/714/e19/tab-figures-data#sec-8
  if (controls) {
    # want to build a separate model for each possible inflection point
    model_stats <- NULL
    for (i in 2:(max(data$months_pre_diag) -1)) {
      data <- data %>% mutate(post_inf = case_when((months_pre_diag <= i) & (case_control == "case") ~ i - months_pre_diag, T ~ as.integer(0)))
      
      # data$post_inf <- 0
      # data$post_inf[(data$months_pre_diag <= i) & (data$case_control == "case")] <- i - data$months_pre_diag[(data$months_pre_diag <= i) & (data$case_control == "case")]
      
      #build model - control for age, sex and year
      # model <- glm.nb(n_consultations ~ months_pre_diag + post_inf + age + female + current_year, data = data)
      model <- glm(n_consultations ~ months_pre_diag + post_inf + age + female + current_year, data = data, family = "poisson")
      #save the model stats
      model_stats <- rbind(model_stats, glance(model))
    }
    
  } else {
    # want to build a separate model for each possible inflection point
    model_stats <- NULL
    
    for (i in 2:(max(data$months_pre_diag) -1)) {
      data <- data %>% mutate(post_inf = case_when((months_pre_diag <= i) ~ i - months_pre_diag, T ~ as.integer(0)))
      # data$post_inf <- 0
      # data$post_inf[data$months_pre_diag <= i] <- i - data$months_pre_diag[data$months_pre_diag <= i] 
      
      #build model - control for age, sex and year
      #model <- glm.nb(n_consultations ~ months_pre_diag + post_inf + age + female + current_year, data = data)
      model <- glm(n_consultations ~ months_pre_diag + post_inf + age + female + current_year, data = data, family = "poisson")
      #save the model stats
      model_stats <- rbind(model_stats, glance(model))
    }
  }
  model_stats$inf_point <- seq(2, max(data$months_pre_diag) - 1)
  
  # find model with largest log-likelihood
  idx <- which.max(model_stats$logLik)
  inflection_point <- model_stats$inf_point[idx]
  
  return(list("inflection_point"=inflection_point))
}

GA_alt_method <- function(data1, controls = F) {
  #alternative - not sure if faster
  #try to vectorise some of the code
  
  #get function to make inf_point columns - only bit which depends on controls
  if (controls) {
    inf_point_function <- function(data2, i) {
      data2 <- data2 %>% mutate("post_inf.{i}" := case_when((months_pre_diag <= i) & (case_control == "case") ~ i - months_pre_diag, T ~ as.integer(0)))
    }
  } else {
    inf_point_function <- function(data2, i) {
      data2 <- data2 %>% mutate("post_inf.{i}" := case_when((months_pre_diag <= i) ~ i - months_pre_diag, T ~ as.integer(0)))
    }
  }
    
  #make each column
  for (i in 2:(max(data1$months_pre_diag) -1)) {
    data1 <- inf_point_function(data1, i)
  }

  #generate each model formula to test - each formula has a different inflection point to test
  model_formulas <- lapply(2:(max(data1$months_pre_diag) -1), function(i) as.formula(sprintf("n_consultations ~ months_pre_diag + post_inf.%d + age + female + current_year", i)))

  #function to build each model and return it's log-likelihood
  model_func <- function(formula, data3 = data1) {
    #model <- glm.nb(formula, data = data3)
    model <- glm(formula, data= data3, family = "poisson")
    logLik(model)
  }
  
  #apply function to every formula
  res <- purrr::map(model_formulas, model_func)

  #find inflection point that maximises logLik
  inflection_point <- (2:(max(data1$months_pre_diag) -1))[which.max(res)]
 
  return(list("inflection_point"=inflection_point)) 
}


#### bootstrap functions ####

#method to use for bootstrapping - based on GA alt and optimised for speed
GA_boot <- function(data1, model_formulas) {
  model_func <- function(formula, data2 = data1) {
    model <- glm(formula, data = data2, family = "poisson")
    logLik(model)
  }
  
  #apply function to every formula
  res <- purrr::map(model_formulas, model_func)
  
  #find inflection point that maximises logLik
  inflection_point <- (2:(max(data1$months_pre_diag) -1))[which.max(res)]
  
  return(list("inflection_point"=inflection_point)) 
}
    
boot_internal_function <- function(data, indices, full_data, model_formulas, controls = F) {
  #here data is the list of patient ids
  #full_data is the full dataset
  if (controls) {
    patient_ids <- data[indices] # allows boot to select data
    bootstrapped_data <- full_data[full_data$matched_case %in% patient_ids,]
    return(GA_boot(bootstrapped_data, model_formulas = model_formulas)$inflection_point)
    
  } else {
    patient_ids <- data[indices] # allows boot to select sample
    bootstrapped_data <- full_data[full_data$patid %in% patient_ids,]
    return(GA_boot(bootstrapped_data, model_formulas = model_formulas)$inflection_point)
  }
}

#this function is set up for parallel processing on a Windows computer - if this is not available some modifications will be needed
bootstrap_DW <- function(data, controls = F, n_reps = 1000) {
  #add inf point columns and generate formulas
  if (controls) {
    inf_point_function <- function(data2, i) {
      data2 <- data2 %>% mutate("post_inf.{i}" := case_when((months_pre_diag <= i) & (case_control == "case") ~ i - months_pre_diag, T ~ as.integer(0)))
    }
  } else {
    inf_point_function <- function(data2, i) {
      data2 <- data2 %>% mutate("post_inf.{i}" := case_when((months_pre_diag <= i) ~ i - months_pre_diag, T ~ as.integer(0)))
    }
  }
  
  #make each column
  for (i in 2:(max(data$months_pre_diag) -1)) {
    data <- inf_point_function(data, i)
  }
  
  #generate each model formula to test - each formula has a different inflection point to test
  model_formulas <- lapply(2:(max(data$months_pre_diag) -1), function(i) as.formula(sprintf("n_consultations ~ months_pre_diag + post_inf.%d + age + female + current_year", i)))
  ncpus <- 8
  cl <- makeCluster(ncpus)
  par.setup <- parLapply(cl, 1:length(ncpus), function(xx){require(stats, dplyr)})
  clusterExport(cl, c('glm', 'logLik', 'GA_boot'))
  
  if (controls) {
    bootstrap <- boot(data = unique(data$patid[data$case_control == "case"]), statistic = boot_internal_function, R = n_reps, full_data = data, model_formulas = model_formulas, controls = T, parallel = "snow", cl = cl, ncpus = ncpus)
  } else {
    bootstrap <- boot(data = unique(data$patid), statistic = boot_internal_function, R = n_reps, full_data = data, model_formulas = model_formulas, controls = F, parallel = "snow", cl = cl, ncpus = ncpus)
  }
  
  stopCluster(cl)
  return(bootstrap)
}

#define the function to apply to each set of parameters
row_function <- function(n_cases, control_ratio, inf_point, max_data_duration, inf_coeff) {
  #set seed
  set.seed(100)
  
  #generate data and apply bootstrap
  if (control_ratio == 0) {
    data <- simulate_data(n_cases = n_cases, inflection_point = inf_point, max_data_duration = max_data_duration, inflection_coefficient = inf_coeff)
    bs_results <- bootstrap_DW(data)
  } else {
    data <- simulate_data(n_cases = n_cases, inflection_point = inf_point, max_data_duration = max_data_duration, controls = T, control_ratio = control_ratio, inflection_coefficient = inf_coeff)
    bs_results <- bootstrap_DW(data, controls = T)
  }
  
  message(paste("t0 =", bs_results$t0))
  return(bs_results$t)
}

vectorized_row_function <- function(x) {
  row_function(x[1], x[2], x[3], x[4], x[5])
}


#### run bootstrap ####

# ideally we would import all parameter values to test and iterate through them
# however, even when optimised this method is slow (proportional to N patients - check speed below)
# the code below can be used to import parameter values and iterate through them but 
# may be preferable to just call specific values

# #read excel file with param values
# params <- xlsx::read.xlsx('N:/Documents/Atlas - synthetic/bootstrap_results/GA.xlsx', sheetIndex = 1)
# 
# bootstrap_results <- NULL
# # do in chunks of five rows so that some results are saved if it breaks halfway through
# its <- ceiling(nrow(params)/5)
# for (i in 1:its) {
#   idx = 1 + (i-1)*5
#   if (i < its) {
#     res <- apply(params[idx:(idx + 4),], 1, vectorized_row_function)
#   } else {
#     res <- apply(params[idx:nrow(params)], 1, vectorized_row_function)
#   }
# 
#   bootstrap_results <- cbind(bootstrap_results, res)
#   message(paste("Done ", i*5))
# }

# if running one at a time, good idea to keep an eye on the time
# system.time({res <- row_function(1000, 0, 8, 12, 0.25)})
# write.table(res, "clipboard", row.names = FALSE, col.names = FALSE)


#### check speed ####

# # how long does it take to do a single iteration?
# #1000 cases, inf_point = 8, max_data_duration = 24
# test_cases <- simulate_data()
# 
# system.time({res.basic <- GA_basic_method(test_cases)}) #1.5 secs
# system.time({res.alt <- GA_alt_method(test_cases)}) # 1.29 secs
# 
# #5 controls per case
# test_controls <- simulate_data(controls = T)
# 
# system.time({res.basic <- GA_basic_method(test_controls)}) #8.75 secs
# system.time({res.alt <- GA_alt_method(test_controls)}) # 8.44 secs


#### examine results ####
bootstrap_results <- xlsx::read.xlsx('N:/Documents/Atlas - synthetic/bootstrap_results/GA.xlsx', sheetIndex = 2) %>%
  rowwise %>%
  mutate(boot_its = list(c_across(starts_with("NA.")))) %>%
  ungroup %>%
  select(c(n_cases, control_ratio, inf_point, max_data_duration, inf_coeff, t0, boot_its))


#because our results are all bootstrapped integer values, we use non-parametric method to estimate CI

non_param_bootstrap_CI <- function(results, t0) {
  tmp <- results - t0
  UL.tmp <- quantile(tmp, c(.975, .025))
  return(t0 - UL.tmp)
}

bootstrap_results$LCI <- 0
bootstrap_results$UCI <- 0

for (i in 1:nrow(bootstrap_results)) {
  if (i %in% c(29, 47, 48)) {
    message(i)
  }
  else{
    tmp_CIs <- non_param_bootstrap_CI(bootstrap_results$boot_its[[i]], bootstrap_results$t0[[i]])
    bootstrap_results$LCI[[i]] <- tmp_CIs[1]
    bootstrap_results$UCI[[i]] <- tmp_CIs[2]
  }
}


