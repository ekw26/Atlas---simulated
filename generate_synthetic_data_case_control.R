# libraries
library(rapportools)
library(MASS)
library(tidyverse)
library(ggplot2)

# functions
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


generate_synthetic_data <- function(n_cases = 1000, inflection_point = 8, max_data_duration = 120, controls = F, control_ratio = 5) {
  # n_cases is number of cases to generate data for
  # inflection_point is number of time units pre-diagnosis that inflection point should occur
  # max_data_duration is max number of time units that a patient should have data for
  # controls is whether to produce "matched" controls
  # control_ratio is the number of controls to produce for each case if controls is True
  
  # first make sure there are no problems with entered values
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
  } else if (!is.boolean(controls)) {
    stop(paste("The selected value of controls (", controls, ") is not possible."))
  } else if (controls) {
    if ((control_ratio%%1 != 0) | (control_ratio <= 0)) {
      stop(paste("Controls have been requested but the selected control_ratio (", control_ratio, ") is not possible."))
    }
  }
  
  # now can start building data
  
  if (controls) {
    #need ids for cases and controls
    n_pats <- (control_ratio + 1)*n_cases
    ids <- seq(1, n_pats)
    
    # patient-level clustering
    # this is random noise included to make the default assumptions 'wrong'
    p_cluster <- rnorm(n_pats, 0, 0.02)
    
    # want to duplicate the existing data (ids and patient-level consistent random effect) so that there are n_months rows per patient
    # also within patient, consider the rows to be 'months' of data
    data <- tibble(patid = rep(ids, each = max_data_duration), case_control = c(rep("case", each = max_data_duration*n_cases), rep("control", each = max_data_duration*n_cases*control_ratio)), p_cluster = rep(p_cluster, each = max_data_duration), months_pre_diag = rep(seq(1, max_data_duration), times = n_pats))
    data$matched_case <- data$patid
    data$matched_case[data$case_control == "control"] <- rep(seq(1, n_cases), each = max_data_duration*control_ratio)
    
    # now need to add a variable to indicate whether before or after inflection_point
    # and how many months post inflection 
    data$alt_month <- 0
    # and how many months post inflection: 1 = inflection_month, n-1 = months past inflection
    data[data$months_pre_diag <= inflection_point, 'alt_month'] <- inflection_point - data[data$months_pre_diag <= inflection_point, 'months_pre_diag'] + 1
    
    # generate a person-month variable for consistency
    data$py <- 1
    
    # simulate some data based on round-ish numbers
    # follow a poisson distribution based around the patient-level mean defined by exp of variables
    data$n_consultations <- 0
    data$n_consultations[data$case_control == "case"] <- rpois(n_cases*max_data_duration, exp(
      -1
      - 0.02*data$months_pre_diag
      + 0.2*data$alt_month
      + rnorm(nrow(data), 0, 0.05)
      + data$p_cluster
      + log(data$py)
    ))
    
    data$n_consultations[data$case_control == "control"] <- rpois(n_cases*control_ratio*max_data_duration, exp(
      -1
      - 0.02*data$months_pre_diag # linear trend for increase over time
      + rnorm(nrow(data), 0, 0.05) # random error for each person/month
      + data$p_cluster # patient-level error
      + log(data$py)
    ))
    
  } else {
    n_pats <- n_cases
    ids <- seq(1, n_pats)
    
    # patient-level clustering
    # this is random noise included to make the default assumptions 'wrong'
    p_cluster <- rnorm(n_pats, 0, 0.02)
    
    # want to duplicate the existing data (ids and patient-level consistent random effect) so that there are n_months rows per patient
    # also within patient, consider the rows to be 'months' of data
    data <- tibble(patid = rep(ids, each = max_data_duration), p_cluster = rep(p_cluster, each = max_data_duration), months_pre_diag = rep(seq(1, max_data_duration), times = n_pats))
    
    # now need to add a variable to indicate whether before or after inflection_point
    # and how many months post inflection 
    data$alt_month <- 0
    # and how many months post inflection: 1 = inflection_month, n-1 = months past inflection
    data[data$months_pre_diag <= inflection_point, 'alt_month'] <- inflection_point - data[data$months_pre_diag <= inflection_point, 'months_pre_diag'] + 1
    
    # generate a person-month variable for consistency
    data$py <- 1
    
    # simulate some data based on round-ish numbers
    # follow a poisson distribution based around the patient-level mean defined by exp of variables
    data$n_consultations <- rpois(nrow(data), exp(
      -1
      - 0.02*data$months_pre_diag
      + 0.2*data$alt_month
      + rnorm(nrow(data), 0, 0.05)
      + data$p_cluster
      + log(data$py)
    ))
  }
  
  return(data)
}


plot_rates <- function(data, controls = F, model = NA) {
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
    summary_data <- data %>% group_by(months_pre_diag) %>% summarise(mean.consults = mean(n_consultations, na.rm = T), sd.consults = sd(n_consultations, na.rm = T), n.consults = n()) %>% mutate(se.consults = sd.consults/sqrt(n.consults), lower.ci = mean.consults - qt(1 - (0.05/2), n.consults - 1)*se.consults, upper.ci = mean.consults + qt(1 - (0.05/2), n.consults - 1)*se.consults)
    ggplot(summary_data, aes(months_pre_diag, mean.consults)) + 
      geom_errorbar(aes(ymin=lower.ci, ymax = upper.ci), width = .1) + 
      geom_line() + geom_point() + 
      geom_smooth(aes(months_pre_diag, n_consultations), data, method = "glm", method.args = list(family = "poisson")) +
      labs(title = "Mean number of consultations per month", x = "Months before diagnosis", y = "Number of consultations") + 
      scale_x_reverse(n.breaks = max(summary_data$months_pre_diag))
  }
}