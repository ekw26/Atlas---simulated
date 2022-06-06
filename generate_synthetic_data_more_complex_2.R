library(wakefield)
library(purrr)
library(ggplot2)

generate_synthetic_more_complex <- function(n_cases = 1000, control_ratio = 5, inflection_point = 24) {
  patients <- tibble(patid = seq(1, n_cases), 
                     pracid = rdunif(n_cases, 100), 
                     case_control = "case", 
                     sex = wakefield::sex(n=n_cases), 
                     age = wakefield::age(n_cases),
                     diag_year = 2019 - rdunif(n_cases, 0, 19), 
                     diag_month = sample(1:12, n_cases, replace = T), 
                     max_follow_up = (diag_year - 2000)*12 + diag_month
                     ) %>%
    mutate(female = case_when(sex == "Female" ~ 1, TRUE ~ 0))
  
  patients$matched_case <- patients$patid
  
  patients <- patients %>%
    slice(rep(1:n(), control_ratio+1)) %>%
    mutate(patid = seq(1, n_cases*(control_ratio + 1)), case_control = c(rep("case", n_cases), rep("control", n_cases*control_ratio)))
  
  #randomly assign each patient follow-up between 1 and their max
  rfu <- function(m) sample(1:m, 1)
  patients <- patients %>% 
    mutate(follow_up = pmap_dbl(list(max_follow_up), rfu)) %>%
    select(-max_follow_up)
  
  #want patient-level effect and practice-level effect
  patients$pat_cluster <- rnorm((control_ratio + 1)*n_cases, 0, 0.02)
  
  prac_effect <- function(i) {
    set.seed(i)
    rnorm(1, 0, 0.005)
  } 
  patients <- patients %>% 
    mutate(prac_cluster = pmap_dbl(list(pracid), prac_effect))
  
  consultations <- tibble(patid = rep(patients$patid, patients$follow_up), pracid = rep(patients$pracid, patients$follow_up), pat_cluster = rep(patients$pat_cluster, patients$follow_up), prac_cluster = rep(patients$prac_cluster, patients$follow_up), female = rep(patients$female, patients$follow_up), age = rep(patients$age, patients$follow_up), case_control = rep(patients$case_control, patients$follow_up), matched_case = rep(patients$matched_case, patients$follow_up), months_pre_diag = sequence(patients$follow_up), diag_year = rep(patients$diag_year, patients$follow_up), diag_month = rep(patients$diag_month, patients$follow_up))
  
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
    mutate(current_year = pmap_dbl(list(diag_month, diag_year, months_pre_diag, "year"), get_current_date)) %>%
    mutate(current_month = pmap_dbl(list(diag_month, diag_year, months_pre_diag, "month"), get_current_date)) %>%
    mutate(alt_month = case_when((months_pre_diag <= inflection_point) & (case_control == "case") ~ inflection_point - months_pre_diag, TRUE ~ 0))
  
  
  # generate a person-month variable for consistency
  consultations$py <- 1
  
  #don't want size of change to be proportional to length of DW
  # i.e. for earlier inflection points need to shrink the alt_month coefficient
  inflection_coefficient <- 2/inflection_point
  
  # simulate some data based on round-ish numbers
  # follow a poisson distribution based around the patient-level mean defined by exp of variables
  consultations$n_consultations <- rpois(nrow(consultations), exp(
    -1
    -1e-3*0.05*consultations$current_year
    - 0.005*consultations$months_pre_diag
    + inflection_coefficient*consultations$alt_month
    + rnorm(nrow(consultations), 0, 0.05)
    + consultations$pat_cluster
    + consultations$prac_cluster
    + 0.2*consultations$female
    + 0.006*consultations$age
    + log(consultations$py)
  ))
  
  return_consultations <- consultations %>% 
    select(patid, female, age, case_control, matched_case, current_month, current_year, months_pre_diag, n_consultations)
  return(return_consultations)
}

summarise_consultations <- function(consultations, max_follow_up = F) {
  #calculate summary data - want mean with CI for each month before diagnosis, for cases and controls
  summary_data <- consultations %>% 
  group_by(case_control, months_pre_diag) %>% 
  summarise(mean.consults = mean(n_consultations, na.rm = T), sd.consults = sd(n_consultations, na.rm = T), n_pat_months = n()) %>% 
  mutate(se.consults = sd.consults/sqrt(n_pat_months), lower.ci = mean.consults - qt(1 - (0.05/2), n_pat_months - 1)*se.consults, upper.ci = mean.consults + qt(1 - (0.05/2), n_pat_months - 1)*se.consults)
  
  if (max_follow_up) {
    summary_data <- summary_data[summary_data$months_pre_diag <= max_follow_up,]
  }
  
  #some months may have zeros/wide CIS because of random follow-up
  ggplot(summary_data, aes(months_pre_diag, mean.consults, col=case_control)) +
  geom_errorbar(aes(ymin=lower.ci, ymax = upper.ci), width = .1) +
  geom_line() + geom_point() +
  labs(title = "Mean number of consultations per month", x = "Months before diagnosis", y = "Number of consultations") +
  scale_x_reverse(n.breaks = max(summary_data$months_pre_diag))
  
}
