#should all throw an error
generate_synthetic_data(max_data_duration = -10)
generate_synthetic_data(max_data_duration = 2.4)
generate_synthetic_data(max_data_duration = 0)

generate_synthetic_data(inflection_point = 7.2)
generate_synthetic_data(inflection_point = -14)

generate_synthetic_data(inflection_point = 8, max_data_duration = -2.7)

#should throw a warning, but continue
generate_synthetic_data(inflection_point = 0)
generate_synthetic_data(inflection_point = 10, max_data_duration = 5)

#should all throw an error
generate_synthetic_data(n_cases = -52)
generate_synthetic_data(n_cases = 751.2)
generate_synthetic_data(controls = "This is a test")
generate_synthetic_data(controls = T, control_ratio = -7.5)


generate_synthetic_data(n_cases = 10, max_data_duration = 10)
generate_synthetic_data(n_cases = 10, max_data_duration = 10, controls = T, control_ratio = 3)
