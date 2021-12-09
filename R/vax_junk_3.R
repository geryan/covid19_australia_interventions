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
    count = sum(count),
    .groups = "drop"
  ) %>% 
  group_by(state, age_class, dose1) %>%
  mutate(
    correction = slider::pslide_dbl(
      .l = list(
        date = date,
        doses = count,
        dose_number = 1
      ),
      .f = immunity_lag_correction,
      .before = Inf
    ) %>%
      unlist,
    marginal_effective_count = correction * count,
    new = slider::slide_dbl(
      .x = count,
      .f = function(.x){
        if(length(.x) == 1){
          .x
        } else {
          .x[2]-.x[1]
        }
      },
      .before = 1
    )
  ) %>%
  arrange(state, age_class, dose1, date) %>%
  ungroup %>%
  mutate(
    marginal_effect = case_when(
      dose1 == "az" ~ efficacy_az_1_dose,
      dose1 == "pf" ~ efficacy_pf_1_dose
    ),
    effective_marginal_efficacy = marginal_effective_count * marginal_effect
  )

cumulative_second_dosed_individuals <- air_collapsed %>%
  filter(!is.na(dose2)) %>%
  group_by(
    state, age_class, date, dose1, dose2
  ) %>%
  summarise(
    count = sum(count),
    .groups = "drop"
  ) %>% 
  group_by(state, age_class, dose1, dose2) %>%
  mutate(
    correction = slider::pslide_dbl(
      .l = list(
        date,
        count,
        2
      ),
      .f = immunity_lag_correction,
      .before = Inf
    ) %>%
      unlist,
    marginal_effective_count = correction * count,
    new = slider::slide_dbl(
      .x = count,
      .f = function(.x){
        if(length(.x) == 1){
          .x
        } else {
          .x[2]-.x[1]
        }
      },
      .before = 1
    )
  ) %>%
  arrange(state, age_class, dose1, dose2, date) %>%
  ungroup %>%
  mutate(
    marginal_effect = case_when(
      dose1 == "az" & dose2 == "az" ~ marginal_az2_azaz,
      dose1 == "az" & dose2 == "pf" ~ marginal_pf2_azpf,
      dose1 == "pf" & dose2 == "pf" ~ marginal_pf2_pfpf,
    ),
    effective_marginal_efficacy = marginal_effective_count * marginal_effect
  )

cumulative_third_dosed_individuals <- air_collapsed %>%
  filter(!is.na(dose3)) %>%
  group_by(
    state, age_class, date, dose1, dose2, dose3
  ) %>%
  summarise(
    count = sum(count),
    .groups = "drop"
  ) %>% 
  group_by(state, age_class, dose1, dose2, dose3) %>%
  mutate(
    correction = slider::pslide_dbl(
      .l = list(
        date,
        count,
        2
      ),
      .f = immunity_lag_correction,
      .before = Inf
    ) %>%
      unlist,
    marginal_effective_count = correction * count,
    new = slider::slide_dbl(
      .x = count,
      .f = function(.x){
        if(length(.x) == 1){
          .x
        } else {
          .x[2]-.x[1]
        }
      },
      .before = 1
    )
  ) %>%
  arrange(state, age_class, dose1, dose2, dose3, date) %>%
  ungroup %>%
  mutate(
    marginal_effect = case_when(
      dose2 == "az" & dose3 == "az" ~ marginal_az3_azazaz,
      dose2 == "az" & dose3 == "pf" ~ marginal_pf3_azazpf,
      dose2 == "pf" & dose3 == "pf" ~ marginal_pf3_pfpfpf,
    ),
    effective_marginal_efficacy = marginal_effective_count * marginal_effect
  )

age_distribution_state <- get_age_distribution_by_state()

cumulative_individuals_any_vaccine <- cumulative_first_dosed_individuals %>%
  group_by(state, age_class, date) %>%
  summarise(
    count = sum(count),
    marginal_effective_count = sum(marginal_effective_count),
    .groups = "drop"
  )

effective_efficacy_data <- bind_rows(
  cumulative_first_dosed_individuals %>%
    select(state, age_class, date, effective_marginal_efficacy),
  cumulative_second_dosed_individuals %>%
    select(state, age_class, date, effective_marginal_efficacy),
  cumulative_third_dosed_individuals %>%
    select(state, age_class, date, effective_marginal_efficacy)
) %>%
  group_by(state, age_class, date) %>%
  summarise(
    effective_total_efficacy = sum(effective_marginal_efficacy),
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
    cumulative_individuals_any_vaccine,
    by = c("state", "age_class", "date")
  )


vaccination_effect <- effective_efficacy_data %>%
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
  rename(
    "effective_any_vaccine" = marginal_effective_count
  ) %>%
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




dose_dates <- unique(vaccination_effect$date)

vaccine_effect_timeseries <- bind_rows(
  vaccination_effect[1,],
  vaccination_effect %>%
    mutate(
      date = date + 6
    )
) %>%
  dplyr::select(
    state,
    date,
    effect = effective_vaccination_transmission_multiplier,
    percent_reduction = effective_vaccination_transmission_reduction_percent
  ) %>%
  full_join(
    y = expand_grid(
      date = seq.Date(
        from = min(dose_dates),
        to = max(dose_dates) + 6,
        by = 1
      ),
      state = states
    ),
    by = c("state", "date")
  ) %>%
  arrange(state, date) %>%
  group_by(state) %>%
  mutate(
    effect = ifelse(
      is.na(effect),
      approx(date, effect, date)$y,
      effect
    ),
    percent_reduction = ifelse(
      is.na(percent_reduction),
      approx(date, percent_reduction, date)$y,
      percent_reduction
    )
  ) %>%
  ungroup




data_date <- max(vaccine_effect_timeseries$date)

saveRDS(
  object = vaccine_effect_timeseries %>%
    dplyr::select(-percent_reduction),
  file = "outputs/vaccine_effect_timeseries.RDS"
)

write_csv(
  vaccine_effect_timeseries,
  file = sprintf(
    "outputs/vaccine_effect_timeseries_%s.csv",
    data_date
  )
)


dpi <- 150
font_size <- 16

ggplot(vaccine_effect_timeseries) +
  geom_line(
    aes(
      x = date,
      y = effect,
      colour = state
    ),
    size = 1.5,
    alpha = 1
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Change in transmission potential",
    col = "State"
  ) +
  scale_x_date(
    breaks = "1 month",
    date_labels = "%b %Y"
  ) +
  ggtitle(
    label = "Vaccination effect",
    subtitle = "Change in transmission potential due to vaccination"
  ) +
  cowplot::theme_cowplot() +
  cowplot::panel_border(remove = TRUE) +
  theme(
    strip.background = element_blank(),
    axis.title.y.right = element_text(vjust = 0.5, angle = 90, size = font_size),
    legend.position = c(0.02, 0.135),
    legend.text = element_text(size = font_size),
    axis.text = element_text(size = font_size),
    plot.title = element_text(size = font_size + 8),
    plot.subtitle = element_text(size = font_size)
  ) +
  scale_colour_manual(
    values = c(
      "darkgray",
      "cornflowerblue",
      "chocolate1",
      "violetred4",
      "red1",
      "darkgreen",
      "darkblue",
      "gold1"
    )
  ) +
  scale_y_continuous(
    position = "right",
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1)
  )

ggsave(
  filename = "outputs/figures/vaccination_effect.png",
  dpi = dpi,
  width = 1500 / dpi,
  height = 1250 / dpi,
  scale = 1.2
)


vaccine_effect_timeseries %>%
  group_by(state) %>%
  mutate(
    delta_week = slider::slide(
      .x = -percent_reduction,
      .f = function(x){
        x[1] - x[7]
      },
      .before = 7
    ) %>%
      unlist
  ) %>%
  ggplot() +
  geom_line(
    aes(
      x = date,
      y = delta_week,
      col = state
    ),
    size = 2
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Change in percentage reduction of transmission potential",
    col = "State"
  ) +
  scale_x_date(
    breaks = "1 month",
    date_labels = "%b %Y"
  ) +
  ggtitle(
    label = "Vaccination effect",
    subtitle = "Change in weekly average percentage reduction in transmission potential due to vaccination"
  ) +
  cowplot::theme_cowplot() +
  cowplot::panel_border(remove = TRUE) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0, face = "bold"),
    axis.title.y.right = element_text(vjust = 0.5, angle = 90),
    panel.spacing = unit(1.2, "lines")
  ) +
  scale_colour_manual(
    values = c(
      "darkgray",
      "cornflowerblue",
      "chocolate1",
      "violetred4",
      "red1",
      "darkgreen",
      "darkblue",
      "gold1"
    )
  )


ggsave(
  filename = "outputs/figures/vaccination_weekly_percent_change_in_tp.png",
  dpi = dpi,
  width = 1500 / dpi,
  height = 1250 / dpi,
  scale = 1.2
)



third_dosed_individuals <- cumulative_third_dosed_individuals %>%
  mutate(
    cumulative_count = count,
    net_count = count,
    marginal_effective_cumulative_count = marginal_effective_count,
    marginal_effective_net_count = marginal_effective_count,
    effective_cumulative_count = case_when(
      dose3 == "pf" & dose2 == "az" ~ prop_az2_azazpf * cumulative_count + prop_pf3_azazpf * marginal_effective_cumulative_count,
      TRUE ~ net_count
    ),
    effective_net_count = case_when(
      dose3 == "pf" & dose2 == "az" ~ prop_az2_azazpf * net_count + prop_pf3_azazpf * marginal_effective_net_count,
      TRUE ~ net_count
    ),
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
    net_count,
    effective_cumulative_count,
    effective_net_count
  )

subtract_third_dosed_individuals <- cumulative_third_dosed_individuals %>%
  group_by(state, age_class, date, dose1, dose2) %>%
  summarise(
    subtract = sum(count),
    .groups = "drop"
  )


second_dosed_individuals <- full_join(
  x = cumulative_second_dosed_individuals,
  y = subtract_third_dosed_individuals,
  by = c("state", "age_class", "date", "dose1", "dose2")
) %>%
  mutate(
    cumulative_count = count,
    net_count = count - subtract,
    marginal_effective_cumulative_count = marginal_effective_count,
    marginal_effective_net_count = marginal_effective_count - subtract,
    effective_cumulative_count = case_when(
      dose1 == "az" & dose2 == "az" ~ prop_az1_azaz * cumulative_count + prop_az2_azaz * marginal_effective_cumulative_count,
      dose1 == "az" & dose2 == "pf" ~ prop_az1_azpf * cumulative_count + prop_pf2_azpf * marginal_effective_cumulative_count,
      dose1 == "pf" & dose2 == "pf" ~ prop_pf1_pfpf * cumulative_count + prop_pf2_pfpf * marginal_effective_cumulative_count
    ),
    effective_net_count = case_when(
      dose1 == "az" & dose2 == "az" ~ prop_az1_azaz * net_count + prop_az2_azaz * marginal_effective_net_count,
      dose1 == "az" & dose2 == "pf" ~ prop_az1_azpf * net_count + prop_pf2_azpf * marginal_effective_net_count,
      dose1 == "pf" & dose2 == "pf" ~ prop_pf1_pfpf * net_count + prop_pf2_pfpf * marginal_effective_net_count
    ),
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
    net_count,
    effective_cumulative_count,
    effective_net_count
  )


subtract_second_dosed_individuals <- cumulative_second_dosed_individuals %>%
  group_by(state, age_class, date, dose1) %>%
  summarise(
    subtract = sum(count),
    .groups = "drop"
  )


first_dosed_individuals <- full_join(
  x = cumulative_first_dosed_individuals,
  y = subtract_second_dosed_individuals,
  by = c("state", "age_class", "date", "dose1")
) %>%
  mutate(
    cumulative_count = count,
    net_count = count - subtract,
    marginal_effective_cumulative_count = marginal_effective_count,
    marginal_effective_net_count = marginal_effective_count - subtract,
    effective_cumulative_count = marginal_effective_cumulative_count,
    effective_net_count = marginal_effective_net_count,
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
    net_count,
    effective_cumulative_count,
    effective_net_count
  )



dose_data_schedule <- bind_rows(
  first_dosed_individuals,
  second_dosed_individuals,
  third_dosed_individuals
)


# check returned net counts are the same as raw net counts
# table should be all TRUE
dose_data_schedule %>%
  full_join(
    air_collapsed %>% rename(raw_net_count = count)
  ) %>%
  mutate(counts_agree = raw_net_count == net_count) %>%
  pull(counts_agree) %>%
  table


dose_data_last_vaccine <- dose_data_schedule %>%
  group_by(state, age_class, date, dose_number, last_vaccine) %>%
  summarise(
    cumulative_count = sum(cumulative_count),
    net_count = sum(net_count),
    effective_cumulative_count = sum(effective_cumulative_count),
    effective_net_count = sum(effective_net_count),
    .groups = "drop"
  ) %>%
  arrange(state, age_class, last_vaccine)

write_csv(
  dose_data_schedule,
  file = sprintf(
    "outputs/effective_dose_data_schedule_%s.csv",
    data_date
  )
)

write_csv(
  dose_data_last_vaccine,
  file = sprintf(
    "outputs/effective_dose_data_last_vaccine_%s.csv",
    data_date
  )
)

# fix this and need to split out moderna from pfizer
# ggplot(
#   effective_dose_data
#   ) +
#   geom_bar(
#     aes(
#       x = date,
#       y = doses,
#       fill = age_class
#     ),
#     stat = "identity"
#   ) +
#   facet_grid(
#     state ~ vaccine + dose_number,
#     scales = "free_y"
#   )


# 
# q_age   <- read_csv(file = "data/vaccinatinon/quantium_vaccination_Rollout/Rollout/dim_age_band.csv")
# q_sa4   <- read_csv(file = "data/vaccinatinon/quantium_vaccination_Rollout/Rollout/dim_sa4.csv")
# q_time  <- read_csv(file = "data/vaccinatinon/quantium_vaccination_Rollout/Rollout/dim_time.csv")
# q_brand <- read_csv(file = "data/vaccinatinon/quantium_vaccination_Rollout/Rollout/dim_vaccine.csv")
# q_vacc  <- read_csv(file = "data/vaccinatinon/quantium_vaccination_Rollout/Rollout/vaccinations.csv")
# 

# 
# effective_dose_data %>%
#   group_by(state, age_class, vaccine, date) %>%
#   summarise(total_vaccinees = sum(doses)) %>%
#   ggplot() +
#   geom_bar(
#     aes(
#       x = date,
#       y = total_vaccinees,
#       fill = age_class
#     ),
#     stat = "identity"
#   ) +
#   facet_grid(
#     state ~ vaccine,
#     scales = "free_y"
#   )
# 
# 
# effective_dose_data %>%
#   filter(state == "ACT", date >= "2021-08-30", vaccine == "pf") %>%
#   ggplot() +
#   geom_bar(
#     aes(
#       x = date,
#       y = doses,
#       fill = age_class
#     ),
#     stat = "identity"
#   ) +
#   facet_grid(
#     age_class ~ dose_number,
#     scales = "free_y"
#   )
# 
# 
# effective_dose_data %>%
#   filter(state == "ACT", date >= "2021-08-30", vaccine == "pf", dose_number == "2") %>%
#   ggplot() +
#   geom_bar(
#     aes(
#       x = date,
#       y = doses,
#       fill = age_class
#     ),
#     stat = "identity"
#   ) +
#   facet_grid(
#     age_class ~ state,
#     scales = "free_y"
#   )
# 
# effective_dose_data %>%
#   filter((state == "ACT" | state == "NT"), date >= "2021-08-30", vaccine == "pf", dose_number == "2") %>%
#   ggplot() +
#   geom_line(
#     aes(
#       x = date,
#       y = doses,
#       col = age_class
#     )
#   ) +
#   facet_wrap(
#     ~state,
#     scales = "free_y"
#   )


