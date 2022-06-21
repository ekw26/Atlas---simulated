#libraries
library(tidyverse)
library(broom)
library(boot)
library(parallel)

# basic DW functions

# sustained non-overlapping CIs
non_overlap_CIs <- function(data, controls = F) {
  # as in https://www.tandfonline.com/doi/full/10.1080/02813432.2022.2057054
  if (controls) {
    # first month CIs for cases do not overlap with CIs of controls
    
    # add a factor variable for months so can model time as continuous and at discrete level each month
    data <- data %>% mutate(month_factor = as.factor(months_pre_diag))
    
    neg_bin <- glm.nb(n_consultations ~ months_pre_diag + case_control + case_control:months_pre_diag + case_control:month_factor, data = data)
    
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
    return(list("inflection_point" = inf_point, "n_consultations_cases" = newdata$n_consultations[newdata$case_control == "case"], "n_consultations_controls" = newdata$n_consultations[newdata$case_control == "control"]))
    
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
    return(list("inflection_point" = inf_point, "n_consultations" = newdata$n_consultations))
  }
}

GA_basic_method <- function(data, controls = F) {
  # as in https://bjgp.org/content/72/714/e19/tab-figures-data#sec-8
  if (controls) {
    # want to build a separate model for each possible inflection point
    # models <- list()
    model_stats <- NULL
    
    #alternative - not sure if faster
    # inf_point_function <- function(data1, i) {
    #   data1 <- data1 %>% mutate("post_inf.{i}" := case_when((months_pre_diag <= i) & (case_control == "case") ~ i - months_pre_diag, T ~ as.integer(0)))
    # }
    # 
    # for (i in 2:(max(tmp_data$months_pre_diag) -1)) {
    #   tmp_data <- inf_point_function(tmp_data, i)
    # }
    # 
    # model_formulas <- lapply(2:(max(tmp_data$months_pre_diag) -1), function(i) as.formula(sprintf("n_consultations ~ months_pre_diag + post_inf.%d + age + female + current_year", i)))
    # 
    # model_func <- function(formula, data = tmp_data) {
    #   model <- glm.nb(formula, data = data)
    #   logLik(model)
    # }
    # 
    # res <- purrr::map(model_formulas, model_func)
    # 
    # inflection_point <- (2:(max(tmp_data$months_pre_diag) -1)[which.max(res)]
    # 
    # for (i in 2:(max(data$months_pre_diag) -1)) {
    #   data$post_inf <- 0
    #   data$post_inf[(data$months_pre_diag <= i) & (data$case_control == "case")] <- i - data$months_pre_diag[(data$months_pre_diag <= i) & (data$case_control == "case")] 
    #   model <- glm.nb(n_consultations ~ months_pre_diag + post_inf + age + female + current_year, data = data)
    #   # models[[i-1]] <- model
    #   model_stats <- rbind(model_stats, glance(model))
    # }
    # model_stats$inf_point <- seq(2, max(data$months_pre_diag) - 1)
    
    # find model with largest log-likelihood
    idx <- which.max(model_stats$logLik)
    inflection_point <- model_stats$inf_point[idx]
    # model <- models[[idx]]
                     
    # estimate consultations from the model
    #newdata <- with(data, expand.grid(months_pre_diag = seq(from = min(data$months_pre_diag), to = max(data$months_pre_diag)), case_control = c("case", "control")))
    #newdata$case_control <- as.factor(newdata$case_control)
    #newdata$post_inf <- 0
    #newdata[(newdata$months_pre_diag <= i) & (newdata$case_control == "case"), 'post_inf'] <- i - newdata[(newdata$months_pre_diag <= i) & (newdata$case_control == "case"), 'months_pre_diag']
    #newdata <- cbind(newdata, predict(model, newdata, type = "link", se.fit = T))
    #newdata <- within(newdata, {
        # n_consultations <- exp(fit)
        # LL <- exp(fit - 1.96 *se.fit)
        # UL <- exp(fit + 1.96*se.fit)})
    
    return(list("inflection_point"=inflection_point))
    
    # return(list("inflection_point" = inflection_point, "model" = model, "n_consultations_cases" = newdata$n_consultations[newdata$case_control == "case"], "n_consultations_controls" = newdata$n_consultations[newdata$case_control == "control"]))
  } else {
    # want to build a separate model for each possible inflection point
    # models <- list()
    model_stats <- NULL
    
    for (i in 2:(max(data$months_pre_diag) -1)) {
      data$post_inf <- 0
      data$post_inf[data$months_pre_diag <= i] <- i - data$months_pre_diag[data$months_pre_diag <= i] 
      model <- glm.nb(n_consultations ~ months_pre_diag + post_inf + age + female + current_year, data = data)
      # models[[i-1]] <- model
      model_stats <- rbind(model_stats, glance(model))
    }
    model_stats$inf_point <- seq(2, max(data$months_pre_diag) - 1)
    
    # find model with largest log-likelihood
    idx <- which.max(model_stats$logLik)
    inflection_point <- model_stats$inf_point[idx]
    # model <- models[[idx]]
     
    # estimate consultations from the model
    #newdata <- with(data, expand.grid(months_pre_diag = seq(from = min(data$months_pre_diag), to = max(data$months_pre_diag))))
    #newdata$post_inf <- 0
    #newdata[(newdata$months_pre_diag <= i), 'post_inf'] <- i - newdata[(newdata$months_pre_diag <= i), 'months_pre_diag']
    #newdata <- cbind(newdata, predict(model, newdata, type = "link", se.fit = T))
    #newdata <- within(newdata, {
      # n_consultations <- exp(fit)
      # LL <- exp(fit - 1.96 *se.fit)
      # UL <- exp(fit + 1.96*se.fit)})
  
    return(list("inflection_point"=inflection_point))
  
  # return(list("inflection_point" = inflection_point, "model" = model, "n_consultations" = newdata$n_consultations))
  }
}
# 
# GA_basic_parallel <- function(data, controls = F) {
#   # as in https://bjgp.org/content/72/714/e19/tab-figures-data#sec-8
#   if (controls) {
#     # want to build a separate model for each possible inflection point
#     # models <- list()
#     model_stats <- NULL
#     cluster <- makeCluster(4)
#     
#     for (i in 2:(max(data$months_pre_diag) -1)) {
#       data$post_inf <- 0
#       data$post_inf[(data$months_pre_diag <= i) & (data$case_control == "case")] <- i - data$months_pre_diag[(data$months_pre_diag <= i) & (data$case_control == "case")] 
#       model <- glm.nb(n_consultations ~ months_pre_diag + post_inf + age + female + current_year, data = data)
#       # models[[i-1]] <- model
#       model_stats <- rbind(model_stats, glance(model))
#     }
#     model_stats$inf_point <- seq(2, max(data$months_pre_diag) - 1)
#     
#     # find model with largest log-likelihood
#     idx <- which.max(model_stats$logLik)
#     inflection_point <- model_stats$inf_point[idx]
#     # model <- models[[idx]]
#     
#     # estimate consultations from the model
#     #newdata <- with(data, expand.grid(months_pre_diag = seq(from = min(data$months_pre_diag), to = max(data$months_pre_diag)), case_control = c("case", "control")))
#     #newdata$case_control <- as.factor(newdata$case_control)
#     #newdata$post_inf <- 0
#     #newdata[(newdata$months_pre_diag <= i) & (newdata$case_control == "case"), 'post_inf'] <- i - newdata[(newdata$months_pre_diag <= i) & (newdata$case_control == "case"), 'months_pre_diag']
#     #newdata <- cbind(newdata, predict(model, newdata, type = "link", se.fit = T))
#     #newdata <- within(newdata, {
#     # n_consultations <- exp(fit)
#     # LL <- exp(fit - 1.96 *se.fit)
#     # UL <- exp(fit + 1.96*se.fit)})
#     
#     return(list("inflection_point"=inflection_point))
#     
#     # return(list("inflection_point" = inflection_point, "model" = model, "n_consultations_cases" = newdata$n_consultations[newdata$case_control == "case"], "n_consultations_controls" = newdata$n_consultations[newdata$case_control == "control"]))
#   } else {
#     # want to build a separate model for each possible inflection point
#     
#     model_boot <- function(i) {
#       data$post_inf <- 0
#       data$post_inf[data$months_pre_diag <= i] <- i - data$months_pre_diag[data$months_pre_diag <= i] 
#       model <- glm.nb(n_consultations ~ months_pre_diag + post_inf + age + female + current_year, data = data)
#       res <- rbind(tibble(), glance(model))
#     }
#     
#     model_stats <- NULL
#     
#     model_stats <- data.table::rbindlist(lapply(seq(2, (max(data$months_pre_diag) -1)), model_boot))
#     
#     model_stats$inf_point <- seq(2, max(data$months_pre_diag) - 1)
#     
#     # find model with largest log-likelihood
#     idx <- which.max(model_stats$logLik)
#     inflection_point <- model_stats$inf_point[idx]
#     
#     return(list("inflection_point"=inflection_point))
#     
#   }
# }

GA_basic_parallel <- function(data, controls = F) {
  # as in https://bjgp.org/content/72/714/e19/tab-figures-data#sec-8
  if (controls) {
    # want to build a separate model for each possible inflection point
    model_boot <- function(i) {
      data$post_inf <- 0
      data$post_inf[(data$months_pre_diag <= i) & (data$case_control == "case")] <- i - data$months_pre_diag[(data$months_pre_diag <= i) & (data$case_control == "case")] 
      model <- glm.nb(n_consultations ~ months_pre_diag + post_inf + age + female + current_year, data = data)
      res <- logLik(model)
    }
    
    model_stats <- NULL
    cl <- makeCluster(4)
    par.setup <- parLapply(cl, 1:length(4), function(xx){ require(stats)})
    clusterExport(cl, c('glm.nb', 'logLik', 'data'))
    
    model_stats <- parLapply(cl, seq(2, (max(data$months_pre_diag) -1)), model_boot)
    
    stopCluster(cl)
    
    # find model with largest log-likelihood
    idx <- which.max(model_stats)
    inflection_point <- seq(2, max(data$months_pre_diag) - 1)[idx]
    
    return(list("inflection_point"=inflection_point))
  
    
  } else {
    # want to build a separate model for each possible inflection point
    
    model_boot <- function(i) {
      data$post_inf <- 0
      data$post_inf[data$months_pre_diag <= i] <- i - data$months_pre_diag[data$months_pre_diag <= i] 
      model <- glm.nb(n_consultations ~ months_pre_diag + post_inf + age + female + current_year, data = data)
      res <- logLik(model)
    }
    
    model_stats <- NULL
    cl <- makeCluster(4)
    par.setup <- parLapply(cl, 1:length(4), function(xx){ require(stats)})
    clusterExport(cl, c('glm.nb', 'logLik', 'data'))
    
    model_stats <- parLapply(cl, seq(2, (max(data$months_pre_diag) -1)), model_boot)
    
    stopCluster(cl)
    
    # find model with largest log-likelihood
    idx <- which.max(model_stats)
    inflection_point <- seq(2, max(data$months_pre_diag) - 1)[idx]
    
    return(list("inflection_point"=inflection_point))
    
  }
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
  if (controls) {
    bootstrap <- boot(data = unique(data$patid[data$case_control == "case"]), statistic = boot_internal_function, R = n_reps, full_data = data, method = method, controls = T, parallel = "multicore", ncpus = 16)
  } else {
    bootstrap <- boot(data = unique(data$patid), statistic = boot_internal_function, R = n_reps, full_data = data, method = method, controls = F, parallel = "multicore", ncpus = 16)
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

