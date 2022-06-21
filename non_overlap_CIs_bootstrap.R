#libraries
library(xlsx)
library(MASS)
library(purrr)
library(tidyverse)
library(boot)
library(ggplot2)

#data simulation code
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
    mutate(follow_up = pmap_dbl(list(max_follow_up), rfu)) %>%
    mutate(follow_up = pmap_dbl(list(follow_up), pairmin)) %>%
    dplyr::select(-max_follow_up)
  
  #add patient-level and practice-level effects
  patients$pat_cluster <- rnorm((control_ratio + 1)*n_cases, 0, patient_var)
  
  prac_effect <- function(i) {
    set.seed(i)
    rnorm(1, 0, prac_var)
  }
  patients <- patients %>%
    mutate(prac_cluster = pmap_dbl(list(pracid), prac_effect),
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
    mutate(current_year = pmap_dbl(list(diag_month, diag_year, months_pre_diag, "year"), get_current_date), 
           current_month = pmap_dbl(list(diag_month, diag_year, months_pre_diag, "month"), get_current_date))
  
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
    mutate(age = pmap_dbl(list(current_month, current_year, diag_year, birth_month, age), get_current_age))
  
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

# function to test
# sustained non-overlapping CIs
non_overlap_CIs <- function(data, controls = F, model = "poisson") {
  # as in https://www.tandfonline.com/doi/full/10.1080/02813432.2022.2057054
  if (controls) {
    # first month CIs for cases do not overlap with CIs of controls
    
    # add a factor variable for months so can model time as continuous and at discrete level each month
    data <- data %>% mutate(month_factor = as.factor(months_pre_diag))
    
    if (model == "negbin") {
      model <- glm.nb(n_consultations ~ months_pre_diag + case_control + case_control:months_pre_diag + case_control:month_factor, data = data)
    } else if (model == "poisson") {
      model <- glm(n_consultations ~ months_pre_diag + case_control + case_control:months_pre_diag + case_control:month_factor, data = data, family = "poisson")
    } else {
      stop("The specified model is not possible - please choose negbin or poisson")
    }
    
    newdata <- with(data, expand.grid(months_pre_diag = seq(from = min(data$months_pre_diag), to = max(data$months_pre_diag)), case_control = c("case", "control")))
    newdata$case_control <- as.factor(newdata$case_control)
    newdata$month_factor <- as.factor(newdata$months_pre_diag)
    newdata <- cbind(newdata, predict(model, newdata, type = "link", se.fit = T))
    newdata <- within(newdata, {
      n_consultations <- exp(fit)
      LL <- exp(fit - 1.96 *se.fit)
      UL <- exp(fit + 1.96*se.fit)
    })
    
    inf_point <- 1
    while (newdata$LL[(newdata$months_pre_diag == inf_point) & (newdata$case_control == "case")] > newdata$UL[(newdata$months_pre_diag == inf_point) & (newdata$case_control == "control")]) {
      inf_point <- inf_point + 1
    }
    
    if (inf_point == max(newdata$months_pre_diag)) {
      inf_point <- 0
      message("No inflection point was detected in the data")
    }
    return(list("inflection_point" = inf_point))
    
  } else {
    # case-only
    # first month CIs do not overlap with CIs of previous month
    
    # fit neg bin regression and get 95% CI for each month
    
    # add a factor variable for months so can model time as continuous and at discrete level each month
    data <- data %>% mutate(month_factor = as.factor(months_pre_diag))
    
    if (model == "negbin") {
      model <- glm.nb(n_consultations ~ months_pre_diag + month_factor, data = data)
    } else if (model == "poisson") {
      model <- glm(n_consultations ~ months_pre_diag + month_factor, data = data, family = "poisson")
    } else {
      stop("The specified model is not possible - please choose negbin or poisson")
    }
    
    newdata <- with(data, expand.grid(months_pre_diag = seq(from = min(data$months_pre_diag), to = max(data$months_pre_diag))))
    newdata$month_factor <- as.factor(newdata$months_pre_diag)
    newdata <- cbind(newdata, predict(model, newdata, type = "link", se.fit = T))
    newdata <- within(newdata, {
      n_consultations <- exp(fit)
      LL <- exp(fit - 1.96 *se.fit)
      UL <- exp(fit + 1.96*se.fit)
    })
    
    inf_point <- 1
    while (newdata$LL[newdata$months_pre_diag == inf_point] > newdata$UL[newdata$months_pre_diag == inf_point + 1]) {
      inf_point <- inf_point + 1
    }
    
    if (inf_point == max(newdata$months_pre_diag)) {
      inf_point <- 0
      message("No inflection point was detected in the data")
    }
    return(list("inflection_point" = inf_point))
  }
}

#bootstrap functions
boot_internal_function <- function(data, indices, full_data, controls = F) {
  #here data is the list of patient ids
  #full_data is the full dataset
  if (controls) {
    patient_ids <- data[indices] # allows boot to select data
    bootstrapped_data <- full_data[full_data$matched_case %in% patient_ids,]
    return(non_overlap_CIs(bootstrapped_data, controls = T)$inflection_point)
    
  } else {
    patient_ids <- data[indices] # allows boot to select sample
    bootstrapped_data <- full_data[full_data$patid %in% patient_ids,]
    return(non_overlap_CIs(bootstrapped_data)$inflection_point)
  }
}

bootstrap_DW <- function(data, controls = F, n_reps = 1000) {
  if (controls) {
    bootstrap <- boot(data = unique(data$patid[data$case_control == "case"]), statistic = boot_internal_function, R = n_reps, full_data = data, controls = T, parallel = "multicore", ncpus = 12)
  } else {
    bootstrap <- boot(data = unique(data$patid), statistic = boot_internal_function, R = n_reps, full_data = data, controls = F, parallel = "multicore", ncpus = 16)
  }
  return(bootstrap)
}

#read excel file with param values
params <- read.xlsx('N:/Documents/Atlas - synthetic/bootstrap_results/non_overlap_cis.xlsx', sheetIndex = 1)

#define the function to apply to each row
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
  
  CIs <- boot.ci(bs_results, type = 'basic')
  #because distribution of bootstrap vals is weird-tailed will use basic/pivotal/empical CI
  # like a percentile distribution but forced to contain t0
  
  inf_point <- bs_results$t0
  LCI <- CIs$basic[4]
  UCI <- CIs$basic[5]
  width <- UCI - LCI
  
  return(list("estimated_inf_point" = inf_point, "LCI" = LCI, "UCI" = UCI, "width" = width))
}

vectorized_row_function <- function(x) {
  row_function(x[1], x[2], x[3], x[4], x[5])
}

its <- ceiling(nrow(params)/5) 
for (i in 1:its) {
  idx = 1 + (i-1)*5
  if (i < its) {
    res <- apply(params[idx:(idx + 4),], 1, vectorized_row_function)
  } else {
    res <- apply(params[idx:nrow(params)], 1, vectorized_row_function)
  }
  
  for (j in 0:4) {
    params[idx + j, "estimated_inf_point"] <- res[[j+1]]$estimated_inf_point
    if (is.null(res[[j+1]]$LCI)) {
      params[idx + j, "LCI"] <- res[[j+1]]$estimated_inf_point
      params[idx + j, "UCI"] <- res[[j+1]]$estimated_inf_point
      params[idx + j, "width"] <- 0
    } else {
      params[idx + j, "LCI"] <- res[[j+1]]$LCI
      params[idx + j, "UCI"] <- res[[j+1]]$UCI
      params[idx + j, "width"] <- res[[j+1]]$width
    }
  }
}



data <- simulate_data(n_cases = 10000)
plots <- c()
for (i in c(50, 100, 500, 1000, 2500, 5000, 10000)) {
  patient_ids <- sample(10000, i, replace = F)
  tmp_data <- data[data$patid %in% patient_ids, ]
  no_inf_point <- non_overlap_CIs(tmp_data)$inflection_point
  ga_inf_point <- GA_basic_parallel(tmp_data)$inflection_point
  message(paste(c('N cases =', i)))
  message(paste(c('NO Inflection point estimate = ', no_inf_point)))
  message(paste(c('GA Inflection point estimate = ', ga_inf_point)))
  
  summary_data <- tmp_data %>% group_by(months_pre_diag) %>% summarise(mean.consults = mean(n_consultations, na.rm = T), sd.consults = sd(n_consultations, na.rm = T), n.consults = n()) %>% mutate(se.consults = sd.consults/sqrt(n.consults), lower.ci = mean.consults - qt(1 - (0.05/2), n.consults - 1)*se.consults, upper.ci = mean.consults + qt(1 - (0.05/2), n.consults - 1)*se.consults)
  g <- ggplot(summary_data, aes(months_pre_diag, mean.consults)) + 
    geom_errorbar(aes(ymin=lower.ci, ymax = upper.ci), width = .1) + 
    geom_line() + geom_point() + 
    geom_vline(xintercept = no_inf_point, colour = "red", linetype = "longdash") +
    geom_vline(xintercept = ga_inf_point, colour = "green", linetype = "longdash") +
    
    geom_smooth(aes(months_pre_diag, n_consultations), data, method = "glm", method.args = list(family = "poisson")) +
    labs(title = paste("Mean number of consultations per month - ", i, "patients"), x = "Months before diagnosis", y = "Number of consultations") + 
    scale_x_reverse(n.breaks = max(summary_data$months_pre_diag))
  plots <- c(plots, list("n_pats" = i, "plot" = g))
}


     