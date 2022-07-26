#libraries
library(tidyverse)
library(broom)
library(boot)
library(parallel)

# basic DW functions

# sustained non-overlapping CIs

non_overlap_CIs <- function(data, controls = F, model = "poisson") {
  # as in https://www.tandfonline.com/doi/full/10.1080/02813432.2022.2057054
  if (controls) {
    # first month CIs for cases do not overlap with CIs of controls
    
    # add a factor variable for months so can model time as continuous and at discrete level each month
    data1 <- data %>% 
      mutate(month_factor = as.factor(months_pre_diag)) %>%
      dplyr::select(month_factor, case_control, n_consultations) %>%
      dplyr::group_by(month_factor, case_control) %>%
      dplyr::summarize(n_consultations = sum(n_consultations), n_patients = dplyr::n()) %>%
      dplyr::ungroup()
    
    if (model == "negbin") {
      model <- glm.nb(n_consultations ~ case_control + case_control:month_factor, data = data1, offset = log(n_patients))
    } else if (model == "poisson") {
      model <- glm(n_consultations ~ case_control + case_control:month_factor, data = data1, family = "poisson", offset= log(n_patients))
    } else {
      stop("The specified model is not possible - please choose negbin or poisson")
    }
    
    # newdata <- with(data, expand.grid(months_pre_diag = seq(from = min(data$months_pre_diag), to = max(data$months_pre_diag)), case_control = c("case", "control")))
    # newdata$case_control <- as.factor(newdata$case_control)
    # newdata$month_factor <- as.factor(newdata$months_pre_diag)
    # newdata <- cbind(newdata, predict(model, newdata, type = "link", se.fit = T))
    # newdata <- within(newdata, {
    #   n_consultations <- exp(fit)
    #   LL <- exp(fit - 1.96*se.fit)
    #   UL <- exp(fit + 1.96*se.fit)
    # })
    # 
    
    newdata <- data1 %>% dplyr::select(-n_consultations)
    newdata <- cbind(newdata, predict(model, newdata, type = "link", se.fit = T)) %>%
      mutate(n_consultations = exp(fit)/n_patients, 
             LL = exp(fit - 1.96*se.fit)/n_patients, 
             UL = exp(fit + 1.96*se.fit)/n_patients)
    
    inf_point <- 1
    while (newdata$LL[(newdata$month_factor == inf_point) & (newdata$case_control == "case")] > newdata$UL[(newdata$month_factor == inf_point) & (newdata$case_control == "control")]) {
      inf_point <- inf_point + 1
    }
    
    if (inf_point == max(data$months_pre_diag)) {
      inf_point <- 0
      message("No inflection point was detected in the data")
    }
    return(list("inflection_point" = inf_point))
    
  } else {
    # case-only
    # first month CIs do not overlap with CIs of previous month
    
    # add a factor variable for months so can model time as continuous and at discrete level each month
    data1 <- data %>% 
      mutate(month_factor = as.factor(months_pre_diag)) %>%
      dplyr::select(month_factor, n_consultations) %>%
      dplyr::group_by(month_factor) %>%
      dplyr::summarize(n_consultations = sum(n_consultations), n_patients = dplyr::n()) %>%
      dplyr::ungroup()
    
    if (model == "negbin") {
      model <- glm.nb(n_consultations ~ month_factor, data = data1, offset = log(n_patients))
    } else if (model == "poisson") {
      model <- glm(n_consultations ~ month_factor, data = data1, family = "poisson", offset = log(n_patients))
    } else {
      stop("The specified model is not possible - please choose negbin or poisson")
    }
    
    # newdata <- with(data, expand.grid(months_pre_diag = seq(from = min(data$months_pre_diag), to = max(data$months_pre_diag))))
    # newdata$month_factor <- as.factor(newdata$months_pre_diag)
    # newdata <- cbind(newdata, predict(model, newdata, type = "link", se.fit = T))
    # newdata <- within(newdata, {
    #   n_consultations <- exp(fit)
    #   LL <- exp(fit - 1.96 *se.fit)
    #   UL <- exp(fit + 1.96*se.fit)
    # })
    newdata <- data1 %>% dplyr::select(-n_consultations)
    newdata <- cbind(newdata, predict(model, newdata, type = "link", se.fit = T)) %>%
      mutate(n_consultations = exp(fit)/n_patients, 
             LL = exp(fit - 1.96*se.fit)/n_patients, 
             UL = exp(fit + 1.96*se.fit)/n_patients)
    
    inf_point <- 1
    while (newdata$LL[newdata$month_factor == inf_point] > newdata$UL[newdata$month_factor == inf_point + 1]) {
      inf_point <- inf_point + 1
    }
    
    if (inf_point == max(data$months_pre_diag)) {
      inf_point <- 0
      message("No inflection point was detected in the data")
    }
    return(list("inflection_point" = inf_point))
  }
}


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

GA_par_method <- function(data, controls = F, ncpus = 4) {
  # as in https://bjgp.org/content/72/714/e19/tab-figures-data#sec-8
  if (controls) {
    # want to build a separate model for each possible inflection point
    
    model_func <- function(i) {
      data <- data %>% mutate(post_inf = case_when((months_pre_diag <= i) & (case_control == "case") ~ i - months_pre_diag, T ~ as.integer(0)))
      # model <- glm.nb(n_consultations ~ months_pre_diag + post_inf + age + female + current_year, data = data)
      model <- glm(n_consultations ~ months_pre_diag + post_inf + age + female + current_year, data = data, family = "poisson")
      res <- logLik(model)
    }
  }
  else {
    
    model_func <- function(i) {
      data <- data %>% mutate(post_inf = case_when((months_pre_diag <= i) ~ i - months_pre_diag, T ~ as.integer(0)))
      #model <- glm.nb(n_consultations ~ months_pre_diag + post_inf + age + female + current_year, data = data)
      model <- glm(n_consultations ~ months_pre_diag + post_inf + age + female + current_year, data = data, family = "poisson")
      res <- logLik(model)
    }
    
  }
  #set up memory for results
  model_stats <- NULL
  
  #create cluster and make sure each core has right packages
  cl <- makeCluster(ncpus)
  par.setup <- parLapply(cl, 1:length(ncpus), function(xx){require(stats, dplyr)})
  clusterExport(cl, c('glm.nb', 'logLik', 'data', '%>%', 'mutate', 'case_when'))
  
  #build all the models
  model_stats <- parLapply(cl, seq(2, (max(data$months_pre_diag) -1)), model_func)
  
  stopCluster(cl)
  
  # find model with largest log-likelihood
  idx <- which.max(model_stats)
  inflection_point <- seq(2, max(data$months_pre_diag) - 1)[idx]
  
  return(list("inflection_point"=inflection_point))
}


GA_agg_method <- function(data1, controls = F) {
  
  # aggregate all data first
  if (controls) {
    data2 <- data1 %>%
      group_by(months_pre_diag, age, female, current_year, case_control) %>%
      summarise(n_consultations = sum(n_consultations), n_patients = n()) %>%
      ungroup()
  } else {
    data2 <- data1 %>%
      group_by(months_pre_diag, age, female, current_year) %>%
      summarise(n_consultations = sum(n_consultations), n_patients = n()) %>% 
      ungroup()
  }
  
  inf_point_function <- function(data3, i) {
    data3 <- data3 %>% mutate("post_inf.{i}" := case_when((months_pre_diag <= i) ~ as.integer(i - months_pre_diag), T ~ as.integer(0)))
  }
  
  #make each column
  for (i in 2:(max(data2$months_pre_diag) -1)) {
    data2 <- inf_point_function(data2, i)
  }
  
  #generate each model formula to test - each formula has a different inflection point to test
  if (controls) {
    model_formulas <- lapply(2:(max(data2$months_pre_diag) -1), function(i) as.formula(sprintf("n_consultations ~ months_pre_diag + case_control:post_inf.%d + age + female + current_year", i)))
  } else {
    model_formulas <- lapply(2:(max(data2$months_pre_diag) -1), function(i) as.formula(sprintf("n_consultations ~ months_pre_diag + post_inf.%d + age + female + current_year", i)))
  }
  
  #function to build each model and return it's log-likelihood
  model_func <- function(formula, data4 = data2) {
    model <- glm(formula, data= data4, family = "poisson", offset = log(n_patients))
    logLik(model)
  }
  
  #apply function to every formula
  res <- purrr::map(model_formulas, model_func)
  
  #find inflection point that maximises logLik
  inflection_point <- (2:(max(data2$months_pre_diag) -1))[which.max(res)]
  
  return(list("inflection_point"=inflection_point)) 
}


summary_stats <- function(vector){
  vec_sd <- sd(vector)
  n <- length(vector)
  vec_mean <- mean(vector)
  error <- qnorm(0.975)*vec_sd/sqrt(n)
  
  return(list("mean" = vec_mean, "UL" = vec_mean + error, "LL" = vec_mean - error))
}

boot_internal_function <- function(data, indices, full_data, method, controls = F) {
  #here data is the list of patient ids
  #full_data is the full dataset
  if (controls) {
    patient_ids <- data[indices] # allows boot to select data
    bootstrapped_data <- full_data[full_data$matched_case %in% patient_ids,]
    return(method(bootstrapped_data, controls = T)$inflection_point)
    
  } else {
    patient_ids <- data[indices] # allows boot to select sample
    bootstrapped_data <- full_data[full_data$patid %in% patient_ids,]
    return(method(bootstrapped_data)$inflection_point)
  }
}

bootstrap_DW_methods <- function(data, controls = F, n_reps = 100, method = non_overlap_CIs) {
  ncpus <- 8
  cl <- makeCluster(ncpus)
  par.setup <- parLapply(cl, 1:length(ncpus), function(xx){require(stats, dplyr)})
  clusterExport(cl, c('glm', 'logLik', 'GA_boot'))
  
  if (controls) {
    bootstrap <- boot(data = unique(data$patid[data$case_control == "case"]), statistic = boot_internal_function, R = n_reps, full_data = data, method = method, controls = T, parallel = "snow", cl = cl, ncpus = ncpus)
  } else {
    bootstrap <- boot(data = unique(data$patid), statistic = boot_internal_function, R = n_reps, full_data = data, method = method, controls = F, parallel = "snow", cl = cl, ncpus = ncpus)
  }
  return(bootstrap)
}



test_DW_methods <- function(n_cases = 1000, inflection_point = 8, max_data_duration = 24, n_reps = 100, method = non_overlap_CIs) {

  tests <- rep(0, n_reps)
  n_consultations <- matrix(0, nrow = n_reps, ncol = max_data_duration)
  for (i in 1:n_reps) {
    results <- method(simulate_data(n_cases = n_cases, inflection_point = inflection_point, max_data_duration = max_data_duration))
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
    results <- non_overlap_CIs(simulate_data(n_cases = n_cases, inflection_point = inflection_point, max_data_duration = max_data_duration, controls = T, control_ratio = control_ratio), controls = T)
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


#test everything

test_cases <- simulate_data(n_cases = 1000, inflection_point = 14, age_range = 30:59)
test_controls <- simulate_data(n_cases = 1000, inflection_point = 4, max_data_duration = 12, controls = T)

bs_case_CI <- bootstrap_DW_methods(test_cases, n_reps = 500)
bs_case_GA <- bootstrap_DW_methods(test_cases, method = GA_basic_method)
bs_case_GA_par <- bootstrap_DW_methods(test_cases, method = GA_basic_parallel)
bs_control_CI <- bootstrap_DW_methods(test_controls, controls = T)
bs_control_GA <- bootstrap_DW_methods(test_controls, controls = T, method = GA_basic_method)


plot(bs_case_CI)
boot.ci(bs_case_CI)

