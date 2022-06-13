#should all throw an error
simulate_data(max_data_duration = -10)
simulate_data(max_data_duration = 2.4)
simulate_data(max_data_duration = 0)

simulate_data(inflection_point = 7.2)
simulate_data(inflection_point = -14)

simulate_data(inflection_point = 8, max_data_duration = -2.7)

#should throw a warning, but continue
simulate_data(inflection_point = 0, n_cases = 5)
simulate_data(inflection_point = 16, max_data_duration = 12, n_cases = 5)

#should all throw an error
simulate_data(n_cases = -52)
simulate_data(n_cases = 751.2)
simulate_data(controls = "This is a test")
simulate_data(controls = T, control_ratio = -7.5)


simulate_data(n_cases = 10, max_data_duration = 10)
simulate_data(n_cases = 10, max_data_duration = 10, controls = T, control_ratio = 3)

test_cases <- simulate_data(n_cases = 1000, inflection_point = 14, age_range = 30:59)
plot_rates(test_cases)

test_controls <- simulate_data(n_cases = 1000, inflection_point = 4, max_data_duration = 12, controls = T)
plot_rates(test_controls, controls = T, model = F)

no_inf_point <- simulate_data(n_cases = 500, inflection_point = 0, controls = T, control_ratio = 3)
plot_rates(no_inf_point)

long <- simulate_data(n_cases = 5000, inflection_point = 60, max_data_duration = 84, controls = T)
plot_rates(long, controls = T, model = F)

before <- simulate_data(n_cases = 500, inflection_point = 36, controls = T)
plot_rates(before, controls = T, model = F)
