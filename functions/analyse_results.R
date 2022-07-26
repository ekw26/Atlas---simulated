library(tidyverse)

#### examine results ####
GA_results <- xlsx::read.xlsx('N:/Documents/Atlas - synthetic/bootstrap_results/GA.xlsx', sheetIndex = 2) %>%
  mutate(controls = (control_ratio > 0), 
         inf_strength = as.factor(inf_coeff*inf_point),
         CI_width = UCI - LCI,
         adjusted_var = var/max_data_duration)

GA_results %>% 
  filter(max_data_duration < 80) %>%
  ggplot(aes(x = n_subjects, y = var)) +
  geom_point(aes(shape = controls, col = inf_strength)) + 
  geom_smooth(aes(col = inf_strength), se = F) +
  facet_wrap( ~ max_data_duration) +
  ggtitle("Variance of bootstrap results against total sample size - Method 2") + 
  xlab("Total sample size") + ylab("Variance")

GA_results %>% 
  filter(max_data_duration < 80) %>%
  ggplot(aes(x = n_cases, y = var)) +
  geom_point(aes(shape = controls, col = inf_strength)) + 
  geom_smooth(aes(col = inf_strength), se = F) +
  facet_grid(controls ~ max_data_duration) +
  ggtitle("Variance of bootstrap results against number of cases - Method 2") + 
  xlab("Number of cases") + ylab("Variance")

GA_results %>% 
  filter(max_data_duration < 80) %>%
  ggplot(aes(x = n_subjects, y = prop_correct)) +
  geom_point(aes(col = inf_strength, size = control_ratio), alpha = 0.5) + 
  geom_smooth(aes(col = inf_strength), se = F) +
  facet_wrap( ~ max_data_duration) +
  ggtitle("Proportion of correct bootstrap estimates against total sample size - Method 2") + 
  xlab("Total sample size") + ylab("Proportion bootstrap estimates correct correct")

GA_results %>% 
  filter(max_data_duration < 80) %>%
  ggplot(aes(x = prop_correct, y = var)) +
  geom_point(aes(col = inf_strength, size = n_subjects), alpha = 0.5) + 
  geom_smooth(aes(col = inf_strength), se = F) +
  facet_wrap( ~ max_data_duration) +
  ggtitle("Proportion of correct bootstrap estimates against variance - Method 2") + 
  ylab("Variance") + xlab("Proportion bootstrap estimates correct")

GA_results %>% 
  filter(max_data_duration < 80) %>%
  ggplot(aes(x = n_cases, y = control_ratio, z = CI_width)) + 
  geom_contour_filled() + 
  facet_grid(inf_strength~max_data_duration)

#### examine results ####
CI_results <- xlsx::read.xlsx('N:/Documents/Atlas - synthetic/bootstrap_results/non_overlap_cis.xlsx', sheetIndex = 2) %>%
  mutate(controls = (control_ratio > 0), 
         inf_strength = as.factor(inf_coeff*inf_point), 
         CI_width = UCI - LCI,
         adjusted_var = var/max_data_duration)
 
CI_results %>% 
  filter(max_data_duration < 80) %>%
  ggplot(aes(x = n_subjects, y = var)) +
  geom_point(aes(shape = controls, col = inf_strength)) + 
  geom_smooth(aes(col = inf_strength), se = F) +
  facet_wrap( ~ max_data_duration) +
  ggtitle("Variance of bootstrap results against total sample size - Method 1") + 
  xlab("Total sample size") + ylab("Variance")

CI_results %>% 
  ggplot(aes(x = n_cases, y = control_ratio, z = CI_width)) + 
  geom_contour_filled() + 
  facet_grid(inf_strength~max_data_duration)


CI_results %>% 
  ggplot(aes(x = prop_correct, y = var)) +
  geom_point(aes(col = inf_strength, size = n_subjects), alpha = 0.5) + 
  geom_smooth(aes(col = inf_strength), se = F) +
  facet_wrap( ~ max_data_duration) +
  ggtitle("Proportion of correct bootstrap estimates against variance - Method 1") + 
  ylab("Variance of bootstrap estimates") + xlab("Proportion bootstrap estimates correct")
