source("R/lib.R")

source("R/functions.R")

air_data <- load_air_data()


air_collapsed <- air_data %>%
  mutate(
    # Here assuming Moderna Spikevax equivalent to Pfizer Comirnaty so collapsing
    across(
      .cols = starts_with("dose"),
      ~ case_when(
        is.na(.) ~ NA_character_,
        . == "COVAST" ~ "az",
        TRUE ~ "pf"
      )
    ),
    # here assuming that in a mixed schedule you continue to get the level of efficacy
    # offered by the more effective vaccine from the point of receiving the more effective one
    dose2 = case_when(
      dose1 == "pf" & dose2 == "az" ~ "pf",
      TRUE ~ dose2
    ),
    dose3 = case_when(
      dose2 == "pf" & dose3 == "az" ~ "pf",
      TRUE ~ dose3
    )
  ) %>%
  group_by(
    state, age_class, date, dose1, dose2, dose3
  ) %>%
  summarise(
    count = sum(count),
    .groups = "drop"
  )

# efficacy of first and second AZ and PF doses
efficacy_az_1_dose <- combine_efficacy(0.46, 0.02)
efficacy_az_2_dose <- combine_efficacy(0.67, 0.36)
efficacy_pf_1_dose <- combine_efficacy(0.57, 0.13)
efficacy_pf_2_dose <- combine_efficacy(0.80, 0.65)

# calculate marginal efficacy in a schedule
# i.e. for a second or third dose, the difference between the last dose and the previous one
marginal_az2_azaz <- efficacy_az_2_dose - efficacy_az_1_dose
marginal_pf2_azpf <- efficacy_pf_2_dose - efficacy_az_1_dose
marginal_pf2_pfpf <- efficacy_pf_2_dose - efficacy_pf_1_dose

marginal_az3_azazaz <- 0
marginal_pf3_azazpf <- efficacy_pf_2_dose - efficacy_az_2_dose
marginal_pf3_pfpfpf <- 0

# calculate the proportion of efficacy at a given dose given by each dose in a schedule
prop_az1_azaz <- efficacy_az_1_dose / efficacy_az_2_dose
prop_az1_azpf <- efficacy_az_1_dose / efficacy_pf_2_dose
prop_pf1_pfpf <- efficacy_pf_1_dose / efficacy_pf_2_dose

prop_az2_azaz <- marginal_az2_azaz / efficacy_az_2_dose
prop_pf2_azpf <- marginal_pf2_azpf / efficacy_pf_2_dose
prop_pf2_pfpf <- marginal_pf2_pfpf / efficacy_pf_2_dose

prop_az2_azazpf <- efficacy_az_2_dose / efficacy_pf_2_dose
prop_pf3_azazpf <- marginal_pf3_azazpf / efficacy_pf_2_dose

cumulative_first_dosed_individuals <- air_collapsed %>%
  group_by(
    state, age_class, date, dose1
  ) %>%
  summarise(
    cumulative_count = sum(count),
    .groups = "drop"
  )

cumulative_second_dosed_individuals <- air_collapsed %>%
  filter(!is.na(dose2)) %>%
  group_by(
    state, age_class, date, dose1, dose2
  ) %>%
  summarise(
    cumulative_count = sum(count),
    .groups = "drop"
  )

cumulative_third_dosed_individuals <- air_collapsed %>%
  filter(!is.na(dose3)) %>%
  group_by(
    state, age_class, date, dose1, dose2, dose3
  ) %>%
  summarise(
    cumulative_count = sum(count),
    .groups = "drop"
  ) 

cumulative_individuals_any_vaccine <- cumulative_first_dosed_individuals %>%
  group_by(state, age_class, date) %>%
  summarise(
    cumulative_count = sum(cumulative_count),
    .groups = "drop"
  )




third_dosed_individuals <- cumulative_third_dosed_individuals %>%
  mutate(
    net_count = cumulative_count,
    dose_number = 3,
    last_vaccine = dose3
  ) %>%
  dplyr::select(
    state,
    age_class,
    date,
    dose_number,
    last_vaccine,
    dose1,
    dose2,
    dose3,
    cumulative_count,
    net_count
  )

subtract_third_dosed_individuals <- cumulative_third_dosed_individuals %>%
  group_by(state, age_class, date, dose1, dose2) %>%
  summarise(
    subtract = sum(cumulative_count),
    .groups = "drop"
  )


second_dosed_individuals <- full_join(
  x = cumulative_second_dosed_individuals,
  y = subtract_third_dosed_individuals,
  by = c("state", "age_class", "date", "dose1", "dose2")
) %>%
  mutate(
    net_count = cumulative_count - subtract,
    dose_number = 2,
    last_vaccine = dose2,
    dose3 = NA
  ) %>%
  dplyr::select(
    state,
    age_class,
    date,
    dose_number,
    last_vaccine,
    dose1,
    dose2,
    dose3,
    cumulative_count,
    net_count
  )


subtract_second_dosed_individuals <- cumulative_second_dosed_individuals %>%
  group_by(state, age_class, date, dose1) %>%
  summarise(
    subtract = sum(cumulative_count),
    .groups = "drop"
  )


first_dosed_individuals <- full_join(
  x = cumulative_first_dosed_individuals,
  y = subtract_second_dosed_individuals,
  by = c("state", "age_class", "date", "dose1")
) %>%
  mutate(
    net_count = cumulative_count - subtract,
    dose_number = 1,
    last_vaccine = dose1,
    dose2 = NA,
    dose3 = NA
  ) %>%
  dplyr::select(
    state,
    age_class,
    date,
    dose_number,
    last_vaccine,
    dose1,
    dose2,
    dose3,
    cumulative_count,
    net_count
  )



dose_data <- bind_rows(
  first_dosed_individuals,
  second_dosed_individuals,
  third_dosed_individuals
)


# check returned net counts are the same as raw net counts
# table should be all TRUE
dose_data %>%
  full_join(
    air_collapsed %>% rename(raw_net_count = count)
  ) %>%
  mutate(counts_agree = raw_net_count == net_count) %>%
  pull(counts_agree) %>%
  table

dose_data %>%
  
  