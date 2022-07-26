# libraries
library(wakefield)
library(purrr)
library(rapportools)
library(MASS)
library(tidyverse)
library(ggplot2)

# functions

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


plot_rates <- function(data, controls = F, model = FALSE) {
  if (controls) {
    #calculate summary data - want mean with CI for each month before diagnosis, for cases and controls
    summary_data <- data %>% group_by(case_control, months_pre_diag) %>% summarise(mean.consults = mean(n_consultations, na.rm = T), sd.consults = sd(n_consultations, na.rm = T), n.consults = n()) %>% mutate(se.consults = sd.consults/sqrt(n.consults), lower.ci = mean.consults - qt(1 - (0.05/2), n.consults - 1)*se.consults, upper.ci = mean.consults + qt(1 - (0.05/2), n.consults - 1)*se.consults)
    
    #fit neg bin to data and get predicted values for each month 
    if (model) {
      newdata <- data.frame(months_pre_diag = rep(seq(from = min(data$months_pre_diag), to = max(data$months_pre_diag)), 2), case_control = rep(c("case", "control"), each = max(data$months_pre_diag)))
      newdata <- cbind(newdata, predict(model, newdata, type = "link", se.fit = T))
      newdata <- within(newdata, {
        n_consultations <- exp(fit)
        LL <- exp(fit - 1.96 *se.fit)
        UL <- exp(fit + 1.96*se.fit)
        })
      
      ggplot(summary_data, aes(months_pre_diag, mean.consults, col=case_control)) + 
        geom_errorbar(aes(ymin=lower.ci, ymax = upper.ci), width = .1) + 
        geom_line() + geom_point() + 
        geom_ribbon(aes(months_pre_diag, n_consultations, ymin = LL, ymax = UL, fill = case_control), newdata, alpha = .25) + 
        geom_line(aes(months_pre_diag, n_consultations, col = case_control), newdata, size = 2) +
        labs(title = "Mean number of consultations per month", x = "Months before diagnosis", y = "Number of consultations") + 
        scale_x_reverse(n.breaks = max(summary_data$months_pre_diag))
    
    } else {
      ggplot(summary_data, aes(months_pre_diag, mean.consults, col=case_control)) +
        geom_errorbar(aes(ymin=lower.ci, ymax = upper.ci), width = .1) +
        geom_line() + geom_point() +
        geom_smooth(aes(months_pre_diag, n_consultations, col = case_control), data, method = "glm", method.args = list(family = "poisson")) +
        labs(title = "Mean number of consultations per month", x = "Months before diagnosis", y = "Number of consultations") +
        scale_x_reverse(n.breaks = max(summary_data$months_pre_diag))
    }
    
    
  } else {
    #calculate summary data - want mean with CI for each month before diagnosis
    summary_data <- data %>% group_by(months_pre_diag) %>% summarise(mean.consults = mean(n_consultations, na.rm = T), sd.consults = sd(n_consultations, na.rm = T), n.patients = n()) %>% mutate(se.consults = sd.consults/sqrt(n.patients), lower.ci = mean.consults - qt(1 - (0.05/2), n.patients - 1)*se.consults, upper.ci = mean.consults + qt(1 - (0.05/2), n.patients - 1)*se.consults)
    ggplot(summary_data, aes(months_pre_diag, mean.consults)) + 
      geom_errorbar(aes(ymin=lower.ci, ymax = upper.ci), width = .1) + 
      geom_line() + geom_point() + 
      geom_smooth(aes(months_pre_diag, n_consultations), data, method = "glm", method.args = list(family = "poisson")) +
      labs(title = "Mean number of consultations per month", x = "Months before diagnosis", y = "Number of consultations") + 
      scale_x_reverse(n.breaks = max(summary_data$months_pre_diag))
  }
}


plot_mean_var <- function(data, controls = F, theta = 1) {
  nb_var <- function(x) x + theta*(x**2)
  if (controls) {
    summary_data <- data %>% group_by(months_pre_diag, case_control) %>% summarise(mean.consults = mean(n_consultations, na.rm = T), sd.consults = sd(n_consultations, na.rm = T), n.patients = n()) %>% mutate(var.consults = sd.consults**2)
    ggplot(summary_data, aes(mean.consults, var.consults, col = case_control)) + 
      geom_point() + 
      geom_abline(slope = 1, intercept = 0, colour = "red") +
      geom_function(fun = nb_var, colour = "blue") +
      scale_x_log10() + 
      scale_y_log10()
  } else {
    summary_data <- data %>% group_by(months_pre_diag) %>% summarise(mean.consults = mean(n_consultations, na.rm = T), sd.consults = sd(n_consultations, na.rm = T), n.patients = n()) %>% mutate(var.consults = sd.consults**2)
    ggplot(summary_data, aes(mean.consults, var.consults)) + 
      geom_point() + 
      geom_abline(slope = 1, intercept = 0, colour = "red") + # poisson slope var = mean
      geom_function(fun = nb_var, colour = "blue") + # neg bin slope var = mean + mean**2
      scale_x_log10() + 
      scale_y_log10()
  }
}
