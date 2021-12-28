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




edd <- dose_data %>%
  arrange(state, age_class, dose1, dose2, dose3) %>%
  group_by(state, age_class, dose1, dose2, dose3) %>%
  mutate(
    new1 = slide_dbl(
      .x = cumulative_count,
      .f = function(x){
       if(length(x) > 1){
         x[2] - x[1]
       } else {
         x[1]
       }
      },
      .before = 1
    ),
    new2 = slide_dbl(
      .x = cumulative_count,
      .f = function(x){
        if(length(x) > 2){
          x[2] - x[1]
        } else if(length(x) > 1){
          x[1]
        } else{
          0
        }
      },
      .before = 2
    ),
    new3 = slide_dbl(
      .x = cumulative_count,
      .f = function(x){
        if(length(x) > 3){
          x[2] - x[1]
        } else if(length(x) > 2) {
          x[1]
        } else {
          0
        }
      },
      .before = 3
    ),
    initial_efficacy = case_when(
      dose_number == 1 ~ 0,
      dose_number == 2 & dose1 == "az" ~ efficacy_az_1_dose,
      dose_number == 2 & dose1 == "pf" ~ efficacy_pf_1_dose,
      dose_number == 3 & dose2 == "az" ~ efficacy_az_2_dose,
      dose_number == 3 & dose2 == "pf" ~ efficacy_pf_2_dose,
    ),
    full_efficacy = case_when(
      dose_number == 1 & dose1 == "az" ~ efficacy_az_1_dose,
      dose_number == 1 & dose1 == "pf" ~ efficacy_pf_1_dose,
      dose_number == 2 & dose2 == "az" ~ efficacy_az_2_dose,
      dose_number == 2 & dose2 == "pf" ~ efficacy_pf_2_dose,
      dose_number == 3 & dose3 == "az" ~ efficacy_az_2_dose,
      dose_number == 3 & dose3 == "pf" ~ efficacy_pf_2_dose,
    ),
    effective_count = case_when(
      dose_number == 1 & new1 + new2 + new3/2 > net_count ~ 0,
      dose_number == 1 ~ net_count - new1 - new2 - new3 + new3 * 0.5,
      TRUE ~ (full_efficacy * (net_count - new1 - new2) + initial_efficacy * (new1 + new2) + new2 * 0.5 *(full_efficacy - initial_efficacy) ) / full_efficacy,
    ),
    total_efficacy = full_efficacy * effective_count
  ) %>%
  ungroup

effective_any_vaccine <- edd %>%
  group_by(state, age_class, date) %>%
  summarise(
    effective_any_vaccine = sum(effective_count)
  )

eedd <- edd %>%
  select(state, age_class, date, effective_count, total_efficacy) %>%
  group_by(state, age_class, date) %>%
  summarise(
    effective_total_efficacy = sum(total_efficacy),
    .groups = "drop"
  ) %>%
  full_join(
    age_distribution_state,
    by = c("state", "age_class")
  ) %>%
  dplyr::select(-fraction) %>%
  mutate(
    effective_average_efficacy_transmission = effective_total_efficacy/pop
  ) %>%
  left_join(
    effective_any_vaccine,
    by = c("state", "age_class", "date")
  )


vaccination_effect_edd <- eedd %>%
  mutate(
    age_group = age_class %>%
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
      ),
    .after = age_class
  ) %>%
  dplyr::select(-age_class) %>%
  mutate(
    effective_coverage_any_vaccine = effective_any_vaccine / pop
  ) %>%
  arrange(
    state, date, age_group
  ) %>%
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


bind_rows(
  vaccination_effect %>% mutate(set = "original"),
  vaccination_effect_edd %>% mutate(set = "effective")
) %>%
  ggplot() +
  geom_line(
    aes(
      x = date,
      y = effective_vaccination_transmission_multiplier,
      col = state,
      alpha = set
    )
  ) +
  scale_alpha_manual(values = c(1, 0.5))

effective_dose_data_write <- edd %>%
  group_by(state, age_class, date, dose_number, last_vaccine) %>%
  summarise(
    cumulative_doses = sum(cumulative_count),
    net_doses = sum(net_count),
    effective_doses = sum(effective_count),
    .groups = "drop"
  ) %>%
  rename("vaccine" = last_vaccine)

write_csv(
  effective_dose_data_write,
  file = sprintf(
    "outputs/effective_dose_data_%s.csv",
    data_date
  )
)

write_csv(effective_dose_data,
          file = "outputs/effective_dose_data_alternate_20211213.csv")
  