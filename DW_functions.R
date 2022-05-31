#libraries
library(tidyverse)
library(broom)

# basic DW functions

# sustained non-overlapping CIs
non_overlap_CIs <- function(data, controls = F) {
  # as in https://www.tandfonline.com/doi/full/10.1080/02813432.2022.2057054
  if (controls) {
    # first month CIs for cases do not overlap with CIs of controls
    
    # add a factor variable for months so can model time as continuous and at discrete level each month
    data <- data %>% mutate(month_factor = as.factor(months_pre_diag))
    
    neg_bin <- glm.nb(n_consultations ~ months_pre_diag + month_factor + case_control + case_control:months_pre_diag + case_control:month_factor, data = data)
    
    newdata <- with(data, expand.grid(months_pre_diag = seq(from = min(data$months_pre_diag), to = max(data$months_pre_diag)), case_control = c("case", "control")))
    newdata$case_control <- as.factor(newdata$case_control)
    newdata$month_factor <- as.factor(newdata$months_pre_diag)
    newdata <- cbind(newdata, predict(neg_bin, newdata, type = "link", se.fit = T))
    newdata <- within(newdata, {
      n_consultations <- exp(fit)
      LL <- exp(fit - 1.96 *se.fit)
      UL <- exp(fit + 1.96*se.fit)
    })
    
    inf_point <- 0
    while (newdata$LL[(newdata$months_pre_diag == inf_point + 1) & (newdata$case_control == "case")] > newdata$UL[(newdata$months_pre_diag == inf_point + 1) & (newdata$case_control == "control")]) {
      inf_point <- inf_point + 1
    }
    
    if (inf_point == max(newdata$months_pre_diag)) {
      inf_point <- 0
      message("No inflection point was detected in the data")
    }
    return(list("inf_point" = inf_point, "n_consultations_cases" = newdata$n_consultations[newdata$case_control == "case"], "n_consultations_controls" = newdata$n_consultations[newdata$case_control == "control"]))
    
  } else {
    # case-only
    # first month CIs do not overlap with CIs of previous month
    
    # fit neg bin regression and get 95% CI for each month
    
    # add a factor variable for months so can model time as continuous and at discrete level each month
    data <- data %>% mutate(month_factor = as.factor(months_pre_diag))
    
    neg_bin <- glm.nb(n_consultations ~ months_pre_diag + month_factor, data = data)
    
    newdata <- with(data, expand.grid(months_pre_diag = seq(from = min(data$months_pre_diag), to = max(data$months_pre_diag))))
    newdata$month_factor <- as.factor(newdata$months_pre_diag)
    newdata <- cbind(newdata, predict(neg_bin, newdata, type = "link", se.fit = T))
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
    return(list("inf_point" = inf_point, "n_consultations" = newdata$n_consultations))
  }
}


GA_method <- function(data, controls = F) {
  # as in https://bjgp.org/content/72/714/e19/tab-figures-data#sec-8
  if (controls) {
    
  } else {
    
    # want to build a separate model for each possible inflection point
    data <- data %>% mutate(month_factor = as.factor(months_pre_diag))
    
    models <- list()
    model_stats <- tibble(null.deviance, df.null, logLik, AIC, BIC, deviance, df.residual, nobs)
    for (i in 2:(max(data$months_pre_diag)-1)) {
      data$post_inf <- 0
      data[data$months_pre_diag <= i, 'post_inf'] <- i - data[data$months_pre_diag <= i, 'months_pre_diag'] + 1
      
      models[[i - 1]] <- glm.nb(n_consultations ~ months_pre_diag + post_inf, data = data)
      model_stats <- bind_rows(model_stats, glance(models[[i-1]]))
    }
    model_stats$inf_point <- seq(2, max(data$months_pre_diag)-1)
    
    # find model with largest log likelihood
    
    
    return(models)
  }
}


summary_stats <- function(vector){
  vec_sd <- sd(vector)
  n <- length(vector)
  vec_mean <- mean(vector)
  error <- qnorm(0.975)*vec_sd/sqrt(n)
  
  return(list("mean" = vec_mean, "UL" = vec_mean + error, "LL" = vec_mean - error))
}


test_DW_methods <- function(n_cases = 1000, inflection_point = 8, max_data_duration = 24, n_reps = 100) {

  tests <- rep(0, n_reps)
  n_consultations <- matrix(0, nrow = n_reps, ncol = max_data_duration)
  for (i in 1:n_reps) {
    results <- non_overlap_CIs(generate_synthetic_data(n_cases = n_cases, inflection_point = inflection_point, max_data_duration = max_data_duration))
    tests[i] <- results$inf_point
    n_consultations[i,] <- results$n_consultations
    
    if (i%%50 == 0) {
      message(paste("Done", i, "iterations"))
    }
    
  }
  
  inflection_points <- summary_stats(tests)
  
  n_consultations <- as_tibble(n_consultations, name_repair = "minimal")
  n_consultations$iteration <- as.factor(seq(1, n_reps))
  n_consultations <- n_consultations %>%
    gather("months_pre_diag", "mean_consults", -iteration) %>%
    mutate(months_pre_diag = as.integer(str_sub(months_pre_diag, 2, -1)))
  
  
  g <- ggplot(n_consultations, aes(months_pre_diag, mean_consults, col = iteration)) + 
    geom_vline(aes(xintercept = inflection_points$mean), col = "red") +
    annotate(geom = "rect", xmin = inflection_points$LL, xmax = inflection_points$UL, ymin = 0, ymax = Inf, fill = alpha("red", 0.3)) + 
    geom_line(alpha = 0.3, show.legend = F) + geom_point(alpha = 0.3, show.legend = F)  +
    labs(title = "Mean number of consultations per month", x = "Months before diagnosis", y = "Number of consultations") + 
    scale_x_reverse(n.breaks = max_data_duration)

  print(g)
  return(inflection_points)
}

test_DW_methods_controls <- function(n_cases = 1000, control_ratio = 5, inflection_point = 8, max_data_duration = 24, n_reps = 100) {
  
  tests <- rep(0, n_reps)
  n_consultations_cases <- matrix(0, nrow = n_reps, ncol = max_data_duration)
  n_consultations_controls <- matrix(0, nrow = n_reps, ncol = max_data_duration)
  for (i in 1:n_reps) {
    results <- non_overlap_CIs(generate_synthetic_data(n_cases = n_cases, inflection_point = inflection_point, max_data_duration = max_data_duration, controls = T, control_ratio = control_ratio), controls = T)
    tests[i] <- results$inf_point
    n_consultations_cases[i,] <- results$n_consultations_cases
    n_consultations_controls[i,] <- results$n_consultations_controls
    
    if (i%%50 == 0) {
      message(paste("Done", i, "iterations"))
    }
    
  }
  
  inflection_points <- summary_stats(tests)
  
  n_consultations_cases <- as_tibble(n_consultations_cases, name_repair = "minimal")
  n_consultations_cases$iteration <- as.factor(seq(1, n_reps))
  n_consultations_cases <- n_consultations_cases %>%
    gather("months_pre_diag", "mean_consults", -iteration) %>%
    mutate(months_pre_diag = as.integer(str_sub(months_pre_diag, 2, -1)))
  n_consultations_cases$case_control <- "case"
  
  n_consultations_controls <- as_tibble(n_consultations_controls, name_repair = "minimal")
  n_consultations_controls$iteration <- as.factor(seq(1, n_reps))
  n_consultations_controls <- n_consultations_controls %>%
    gather("months_pre_diag", "mean_consults", -iteration) %>%
    mutate(months_pre_diag = as.integer(str_sub(months_pre_diag, 2, -1)))
  n_consultations_controls$case_control <- "control"
  
  n_consultations <- bind_rows(n_consultations_cases, n_consultations_controls)
  
  g <- ggplot(n_consultations, aes(months_pre_diag, mean_consults, group = iteration)) + 
    geom_vline(aes(xintercept = inflection_points$mean), col = "red") +
    annotate(geom = "rect", xmin = inflection_points$LL, xmax = inflection_points$UL, ymin = 0, ymax = Inf, fill = alpha("red", 0.3)) + 
    geom_line(aes(linetype = case_control, col = case_control), alpha = 0.3, show.legend = F) + geom_point(alpha = 0.3, show.legend = F)  +
    labs(title = "Mean number of consultations per month", x = "Months before diagnosis", y = "Number of consultations") + 
    scale_x_reverse(n.breaks = max_data_duration)
  
  print(g)
  return(inflection_points)
}
