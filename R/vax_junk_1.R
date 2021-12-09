# effective dose working
tibble(
  date = rep(
    seq.Date(
      from = as.Date("2021-02-22"),
      to = as.Date("2021-03-29"),
      by = 7
    ),
    times = 2
  ),
  doses = c(10, 20, 20, 20, 10,  0,
            0,  0, 10, 20, 30, 40),
  dose_number = c(
    rep(1, 6),
    rep(2, 6)
  )
) %>%
  group_by(dose_number) %>%
  mutate(
    correction = slider::pslide_dbl(
      .l = list(
        date,
        doses,
        dose_number
      ),
      .f = immunity_lag_correction,
      .before = Inf
    ) %>%
      unlist,
    marginal_effective_doses = correction * doses
  ) %>%
  mutate(
    effective_doses = 0.7*doses + 0.3*marginal_effective_doses
  )

efficacy_az_1_dose <- combine_efficacy(0.46, 0.02)
efficacy_az_2_dose <- combine_efficacy(0.67, 0.36)
efficacy_pf_1_dose <- combine_efficacy(0.57, 0.13)
efficacy_pf_2_dose <- combine_efficacy(0.80, 0.65)

marginal_az_az <- efficacy_az_2_dose - efficacy_az_1_dose
marginal_az_pf <- efficacy_pf_2_dose - efficacy_az_1_dose
marginal_pf_pf <- efficacy_pf_2_dose - efficacy_pf_1_dose

marginal_az_az_az <- 0
marginal_az_az_pf <- efficacy_pf_2_dose - efficacy_az_2_dose
marginal_pf_pf_pf <- 0

# calculate the proportion of efficacy at a given dose
# given by each dose
prop_az1_azaz <- efficacy_az_1_dose / efficacy_az_2_dose
prop_az1_azpf <- efficacy_az_1_dose / efficacy_pf_2_dose
prop_pf1_pfpf <- efficacy_pf_1_dose / efficacy_pf_2_dose

prop_az2_azaz <- marginal_az_az / efficacy_az_2_dose
prop_pf2_azpf <- marginal_az_pf / efficacy_pf_2_dose
prop_pf2_pfpf <- marginal_pf_pf / efficacy_pf_2_dose

prop_az2_azazpf <- efficacy_az_2_dose / efficacy_pf_2_dose
prop_pf3_azazpf <- marginal_az_az_pf / efficacy_pf_2_dose


net_third_dosed_individuals <- cumulative_third_dosed_individuals %>%
  dplyr::select(state, age_class, date, dose1, dose2, dose3, count, effective_count) %>%
  rename(
    net_count = count,
    marginal_effective_net_count = effective_count
  ) %>%
  mutate(
    effective_net_count = case_when(
      dose3 == "pf" & dose2 == "az" ~ prop_az2_azazpf * net_count + prop_pf3_azazpf * marginal_effective_net_count,
      TRUE ~ net_count
    )
  )

subtract_third_dosed_individuals <- cumulative_third_dosed_individuals %>%
  group_by(state, age_class, date, dose1, dose2) %>%
  summarise(
    subtract = sum(count),
    .groups = "drop"
  )


net_second_dosed_individuals <- full_join(
  x = cumulative_second_dosed_individuals,
  y = subtract_third_dosed_individuals,
  by = c("state", "age_class", "date", "dose1", "dose2")
  ) %>%
  mutate(
    net_count = count - subtract,
    marginal_effective_net_count = effective_count - subtract,
    effective_net_count = case_when(
      dose1 == "az" & dose2 == "az" ~ prop_az1_azaz * net_count + prop_az2_azaz * marginal_effective_net_count,
      dose1 == "az" & dose2 == "pf" ~ prop_az1_azpf * net_count + prop_pf2_azpf * marginal_effective_net_count,
      dose1 == "pf" & dose2 == "pf" ~ prop_pf1_pfpf * net_count + prop_pf2_pfpf * marginal_effective_net_count
    )
  ) %>%
  dplyr::select(state, age_class, date, dose1, dose2, net_count, marginal_effective_net_count, effective_net_count)


subtract_second_dosed_individuals <- cumulative_second_dosed_individuals %>%
  group_by(state, age_class, date, dose1) %>%
  summarise(
    subtract = sum(count),
    .groups = "drop"
  )


net_first_dosed_individuals <- full_join(
  x = cumulative_first_dosed_individuals,
  y = subtract_second_dosed_individuals,
  by = c("state", "age_class", "date", "dose1")
) %>%
  mutate(
    net_count = count - subtract,
    marginal_effective_net_count = effective_count - subtract,
    effective_net_count = marginal_effective_net_count
  ) %>%
  dplyr::select(state, age_class, date, dose1, net_count, marginal_effective_net_count, effective_net_count)

################



dd <- dose_data_last_vaccine %>%
  mutate(
    dose_number = ifelse(dose_number == 1, 1, 2)
  ) %>%
  group_by(state, age_class, date, dose_number, last_vaccine) %>%
  summarise(
    doses = sum(net_count),
    effective_doses = sum(effective_net_count),
    .groups = "drop"
  ) %>%
  full_join(
    age_distribution_state,
    by = c("state", "age_class")
  ) %>%
  group_by(state, age_class, date) %>%
  mutate(
    any_vaccine = sum(doses),
    effective_any_vaccine = sum(effective_doses),
  ) %>%
  ungroup %>%
  mutate(
    fraction = doses / any_vaccine,
    effective_fraction = effective_doses / effective_any_vaccine,
    coverage_any_vaccine = any_vaccine / pop,
    effective_coverage_any_vaccine = effective_any_vaccine / pop,
    fraction = ifelse(is.nan(fraction), 0, fraction),
    effective_fraction = ifelse(is.nan(effective_fraction), 0, effective_fraction)
  ) %>%
  arrange(state, age_class, date)

efficacy_data_alt <- dd %>%
  pivot_wider(
    names_from = c(last_vaccine, dose_number),
    values_from = c(doses, effective_doses, fraction, effective_fraction)
  ) %>%
  mutate(
    effective_average_efficacy_transmission = average_efficacy(
      efficacy_az_1_dose = combine_efficacy(0.46, 0.02),
      efficacy_az_2_dose = combine_efficacy(0.67, 0.36),
      efficacy_pf_1_dose = combine_efficacy(0.57, 0.13),
      efficacy_pf_2_dose = combine_efficacy(0.80, 0.65),
      proportion_pf_2_dose = effective_fraction_pf_2,
      proportion_az_2_dose = effective_fraction_az_2,
      proportion_pf_1_dose = effective_fraction_pf_1,
      proportion_az_1_dose = effective_fraction_az_1
    ),
    effective_average_efficacy_transmission = replace_na(effective_average_efficacy_transmission, 0)
  ) %>%
  left_join(
    age_lookup,
    by = c("age_class" = "age_5y")
  )  %>%
  dplyr::select(
    -age
  ) %>%
  rename(
    age = age_class
  )

efficacy_data_alt %>% glimpse

vaccination_effect_alt <- efficacy_data_alt %>%
  mutate(
    age_group = age %>%
      factor(
        levels = c(
          "0-4",
          "5-9",
          "10-14",
          "15-19",
          "20-24",
          "25-29",
          "30-34",
          "35-39",
          "40-44",
          "45-49",
          "50-54",
          "55-59",
          "60-64",
          "65-69",
          "70-74",
          "75-79",
          "80+"
        )
      )
  ) %>%
  arrange(
    state, date, age_group
  )  %>%
  group_by(
    state, date
  ) %>%
  summarise(
    effective_vaccination_transmission_multiplier = vaccination_transmission_effect(
      age_coverage = effective_coverage_any_vaccine,
      efficacy_mean = effective_average_efficacy_transmission,
      next_generation_matrix = baseline_matrix()
    )$overall,
    .groups = "drop"
  ) %>%
  mutate(
    effective_vaccination_transmission_reduction_percent =
      100 * (1 - effective_vaccination_transmission_multiplier)
  ) %>%
  mutate(
    dubious = (date - min(date)) < 21,
    across(
      starts_with("effective_"),
      ~ ifelse(dubious, NA, .)
    )
  ) %>%
  select(
    -dubious
  ) %>% 
  mutate(
    effective_vaccination_transmission_multiplier = ifelse(
      is.na(effective_vaccination_transmission_multiplier),
      1,
      effective_vaccination_transmission_multiplier
    ),
    effective_vaccination_transmission_reduction_percent = ifelse(
      is.na(effective_vaccination_transmission_reduction_percent),
      0,
      effective_vaccination_transmission_reduction_percent
    )
  )


##########


dose_data_schedule %>%
  select(-effective_cumulative_count, -effective_net_count) %>%
  mutate(
    schedule = paste0(
      dose1,
      ifelse(is.na(dose2), "", dose2),
      ifelse(is.na(dose3), "", dose3)
    )
  )


